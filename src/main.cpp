
#if defined(_MSC_VER)
#define NOMINMAX
#endif

#include <iostream>
#include <thread>
#include <chrono>
#include <atomic>

#include <cstddef>
#include <cstdlib>
#include <ctgmath>

#include <omp.h>

#include "SRTM.h"
#include "Mesh.h"
#include "Parsing.h"
#include "Utilities.h"
#include "CLI.h"

template<class T>
Mesh meshFromSrtmTile(const T &tile, const LCS &coordinateSystem)
{
    Mesh mesh;

    auto &coordinates = tile.coordinates();
    mesh.vertexPositions.resize(tile.numValues());
    std::transform(coordinates.begin(), coordinates.end(), mesh.vertexPositions.begin(),
                   [&](WGS84::Point coord) -> Eigen::Vector3f {
                       return coordinateSystem.fromECEF(coord).cast<float>();
                   });

    auto size = tile.size();
    mesh.faces.clear();
    mesh.faces.reserve(2*(size[0]-1)*(size[1]-1));
    for (int x = 0; x < size[0]-1; ++x)
    {
        for (int y = 0; y < size[1]-1; ++y)
        {
            auto idx1 = tile.indexAtPosition({x  ,y  });
            auto idx2 = tile.indexAtPosition({x+1,y  });
            auto idx3 = tile.indexAtPosition({x  ,y+1});
            auto idx4 = tile.indexAtPosition({x+1,y+1});

            mesh.faces.push_back({idx1,idx2,idx3});
            mesh.faces.push_back({idx2,idx4,idx3});
        }
    }

    return mesh;
}

template<class T>
void generateMeshes(const CommandLineArguments &arguments)
{
    typedef T Tile;
    typedef typename Tile::Size Size;
    typedef typename Tile::Bounds Bounds;
    typedef typename Tile::AbstractCache AbstractCache;
    typedef typename Tile::Cache Cache;
    typedef typename Tile::Cache WeakCache;
	
	auto &from   = arguments.boundaryPoints[0];
	auto &to     = arguments.boundaryPoints[1];
	auto &origin = arguments.meshOrigin;

	auto mapBoundaries = Tile::bounds(from,to);
    Eigen::Vector2i mapSize = mapBoundaries.diagonal() + Eigen::Vector2i(1,1);

	Eigen::Vector2i tileSize = arguments.tileSize;
	if ((tileSize.array() < 1).any())
		tileSize = mapSize;

    omp_set_nested(0);
	omp_set_num_threads(arguments.numThreads);

    static const auto ceil = [](double a){ return std::ceil(a); };
    Eigen::Vector2i numTiles = Eigen::Vector2d( //Element-wise division -> ceiling
                                    (mapSize.cast<double>().array()/tileSize.cast<double>().array()
                               ).unaryExpr(ceil)).cast<int>();

	
	std::cout << "Input  directory: " << arguments.inputDirectory << "\n";
    std::cout << "Output directory: " << arguments.outputDirectory << "\n\n";

    std::cout << "From:      " << std::setw(5) << from[0]      << " / " << std::setw(5) << from[1]   << " / " << std::setw(5) << from[2]   << "\n";
    std::cout << "To:        " << std::setw(5) << to[0]        << " / " << std::setw(5) << to[1]     << " / " << std::setw(5) << to[2]     << "\n";
    std::cout << "Origin:    " << std::setw(5) << origin[0]    << " / " << std::setw(5) << origin[1] << " / " << std::setw(5) << origin[2] << "\n\n";

    std::cout << "Tile Size: " << std::setw(5) << tileSize[0]  << " / " << std::setw(5) << tileSize[1] << "\n";
    std::cout << "#Tiles:    " << std::setw(5) << numTiles[0]  << " / " << std::setw(5) << numTiles[1] << "\n\n";

    Tensor<Bounds, 2> tileDefinitions(numTiles);
    for (int y = 0; y < numTiles[1]; ++y)
    {
        for (int x = 0; x < numTiles[0]; ++x)
        {
            Eigen::Vector2i position = (tileSize.array() -1) * Eigen::Array2i(x,y);
            Size size = tileSize.array().min(Eigen::Array2i(mapSize - position));
            tileDefinitions[{x,y}] = Tile::bounds(mapBoundaries.min()+position, size);
        }
    }

    LCS lcs(origin);
//    std::unique_ptr<AbstractCache> cache(new Cache(inputDir));
	std::unique_ptr<AbstractCache> cache(new WeakCache(arguments.inputDirectory));

    const size_t numTasks = tileDefinitions.numValues();
    std::atomic<size_t> tasksCompleted;
    auto start = std::chrono::high_resolution_clock::now();
    tasksCompleted.store(0);

    static const auto print_progress = [&]() -> void {
        size_t lastCounter    = 0;
        size_t currentCounter = 0;
        do
        {
            currentCounter = tasksCompleted.load();
            if (currentCounter != lastCounter)
            {
                auto now = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::seconds>(now-start).count();

                lastCounter = currentCounter;
                std::cout << "Progress: " << std::setw(7) << std::fixed << std::setprecision(3)
                          << ((100.0 * currentCounter) / numTasks) << "% "
                          << "( " << std::setw(5) << currentCounter << " / " << std::setw(5) << numTasks << " completed; "
                          << (((double)currentCounter)/duration) << " Tiles/s)\n";
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
        } while(currentCounter != numTasks);
    };

    static const auto generate_tiles = [&]() -> void {
        //Parallelize the computation
        //Scheduling static,1 will interleave the different threads
        //(i.e. thread 1 gets tile 1, thread 2 tile 2 etc.)
        //This should lead them to sharing more tiles, which means
        //less I/O overhead while reading the tiles from disk.
        #pragma omp parallel for default(shared)
        for (int i = 0; i < tileDefinitions.numValues(); ++i)
        {
            try
            {
                std::stringstream ss;
				ss << arguments.outputDirectory << "/" << i << ".ply";

                auto filename = ss.str();
                //ScopedTimer timer(filename);

                auto &bounds = tileDefinitions[i];
                auto tile = Tile::stitch(bounds,*cache);
                auto mesh = meshFromSrtmTile(tile, lcs);

                std::ofstream file(filename);

				if (arguments.produceAsciiPly)
					mesh.toAsciiPly(file);
				else
					mesh.toBinaryPly(file);
            }
            catch(std::exception &e)
            { //Exceptions break the omp loop and lead to a terminate(). Handle them here.
                std::cerr << "An unexpected error occurred: \n"
                          << e.what() << std::endl;
                exit(-1);
            }
            catch(...)
            {
                std::cerr << "An unexpected error occurred (unknown type)." << std::endl;
                exit(-1);
            }
            ++tasksCompleted;
        }
    };

    auto progressThread   = std::thread(print_progress);
    auto generationThread = std::thread(generate_tiles);

    progressThread.join();
    generationThread.join();

    auto now = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(now-start).count();
    std::cout << "\n"
              << "Elapsed time: " << duration << "s."
              << std::endl;
}


int main(int argc, const char** argv)
{
    try
    {
		auto arguments = parseCommandLine(argc, argv);

		if (arguments.useCoarseSrtmResolution)
			generateMeshes<SRTM3::Tile>(arguments);
		else
			generateMeshes<SRTM1::Tile>(arguments);
    } 
	catch(std::exception &e)
    {
        std::cerr << "An unexpected error occurred: \n"
                  << e.what() << std::endl;
        exit(-1);
    }
}
