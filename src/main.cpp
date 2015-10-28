
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

#include <dirent.h>
#ifdef _MSC_VER
#include <direct.h>
#define getcwd _getcwd
#else
#include <unistd.h>
#endif // _MSC_VER

#include <omp.h>

#include <docopt.h>

#include "SRTM.h"
#include "Mesh.h"
#include "Parsing.h"
#include "Utilities.h"

static const char USAGE[] =
R"(
Create meshes from SRTM files.

Usage:
   srtm2ply [options] <p1> <p2> [<INPUT_DIRECTORY> [<OUTPUT_DIRECTORY>]]

Options:
   -h --help                         Print this page.
   -v --version                      Display the application's version number.
   -o --origin=<origin>              Define the origin of the generated tiles (in WGS84 coordinates).
   -s --tile-size=wxh                Set the maximum size of tiles (in sampling points).
   -a --ascii-ply                    Generate ASCII-PLY meshes.
   -j <n>                            Generate N meshes in parallel.  [default: #CPUs]
   --srtm3                           The input directory contains SRTM data using a 90m resolution.

p1 and p2 denote two diagonally located corners of the bounding box of the desired mesh (using the WGS84
coordinate system), and have to be specified as WGS84 coordinates (as well as the desired origin).

The generated PLY-File will contain a mesh consisting only out of triangles. The data type for the vertex
indices in the PLY-File will be chosen according to the maximum amount of possible points in the mesh.
Thus, the -s option should be chosen accordingly to generate indices of a certain size (i.e. a tile size
of 256x256 sampling points will result in 16 bit indices).

WGS84 coordinates have to be specified using the following format (using degrees as the unit of angle):
<latitude>[N|S] [-]<longitude>[E|W] [<altitude>m], e.g. "51.48125N 0.008069W" for Greenwich.

Tile sizes should be specified as <width>x<height> in integer numbers, with 1 representing the distance between
two sampling points in the SRTM dataset (i.e. an arc second for the newest data).
)";

void printArguments(const std::map<std::string, docopt::value> &args)
{
    std::cout << "Arguments:\n\n";
    for (auto &entry: args)
    {
        const std::string &name = entry.first;
        const docopt::value &value = entry.second;

        std::cout << name << ": \n";

        if (value.isString())
            std::cout << "       (string)   " << value.asString() << "\n\n";
        else if (value.isLong())
            std::cout << "       (long)     " << value.asLong() << "\n\n";
        else if (value.isBool())
            std::cout << "       (bool)     " << value.asBool() << "\n\n";
        else if (value.isStringList())
        {
            const auto &list = value.asStringList();
            std::cout << "       (str_list) [\n";

            for (auto &str: list)
            std::cout << "                        " << str << ",\n";
            std::cout << "                  ]\n\n";
        }
        else if (!value) std::cout << "       (empty)\n\n";
        else std::cout << "       (unknown)\n\n";

    }
    std::cout << std::endl;
}

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
void generateMeshes(std::map<std::string, docopt::value> &args)
{
    typedef T Tile;
    typedef typename Tile::Size Size;
    typedef typename Tile::Bounds Bounds;
    typedef typename Tile::AbstractCache AbstractCache;
    typedef typename Tile::Cache Cache;
    typedef typename Tile::Cache WeakCache;

    const std::string cwd(getcwd(nullptr,0));

    std::string inputDir  = args["<INPUT_DIRECTORY>" ].isString() ? args["<INPUT_DIRECTORY>" ].asString() : cwd;
    std::string outputDir = args["<OUTPUT_DIRECTORY>"].isString() ? args["<OUTPUT_DIRECTORY>"].asString() : cwd;

    inputDir  = absolutePath(inputDir.c_str());
    outputDir = absolutePath(outputDir.c_str());

    WGS84::Point from, to, origin;
    from = parseWGS84(args["<p1>"].asString());
    to   = parseWGS84(args["<p2>"].asString());

    if (args["<origin>"].isString())
        origin = parseWGS84(args["<origin>"].asString());
    else
        origin = WGS84::Point(Eigen::Vector3d((from+to)/2.0));

    auto mapBoundaries = Tile::bounds(from,to);
    Eigen::Vector2i mapSize = mapBoundaries.diagonal() + Eigen::Vector2i(1,1);

    Eigen::Vector2i tileSize;
    if (args["--tile-size"].isString())
        tileSize = parseSize(args["--tile-size"].asString());
    else
        tileSize = mapSize;

    omp_set_nested(0);
    {
        int num_threads = parseNumberOfThreads(args["-j"].asString());
        omp_set_num_threads(num_threads);
    }

    static const auto ceil = [](double a){ return std::ceil(a); };
    Eigen::Vector2i numTiles = Eigen::Vector2d( //Element-wise division -> ceiling
                                    (mapSize.cast<double>().array()/tileSize.cast<double>().array()
                               ).unaryExpr(ceil)).cast<int>();

    std::cout << "Input  directory: " << inputDir  << "\n";
    std::cout << "Output directory: " << outputDir << "\n\n";

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
    std::unique_ptr<AbstractCache> cache(new WeakCache(inputDir));

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
                ss << outputDir << "/" << i << ".ply";

                auto filename = ss.str();
                //ScopedTimer timer(filename);

                auto &bounds = tileDefinitions[i];
                auto tile = Tile::stitch(bounds,*cache);
                auto mesh = meshFromSrtmTile(tile, lcs);

                std::ofstream file(filename);

                if (args["--ascii-ply"] && args["--ascii-ply"].asBool())
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

    auto progressThread = std::thread(print_progress);
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

        auto args = docopt::docopt(USAGE,
                                   { argv + 1, argv + argc },
                                   true,                       // show help if requested
                                   "SRTM2PLY 0.1"              // version string
                                   );

        //printArguments(args);

        if (args["--srtm3"] && args["--srtm3"].asBool())
            generateMeshes<SRTM3::Tile>(args);
        else
            generateMeshes<SRTM1::Tile>(args);

    } catch(std::exception &e)
    {
        std::cerr << "An unexpected error occurred: \n"
                  << e.what() << std::endl;
        exit(-1);
    }
}
