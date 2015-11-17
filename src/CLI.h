
#include <array>
#include <Eigen/Geometry>

#include "Geometry.h"

static const char VERSION[] = "SRTM2PLY 0.1";

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

struct CommandLineArguments
{
	CommandLineArguments():
		inputDirectory("."),
		outputDirectory("."),

		tileSize(0,0),
		useCoarseSrtmResolution(false),
		produceAsciiPly(false),
		numThreads(0)
	{}

	std::string inputDirectory;
	std::string outputDirectory;
	std::array<WGS84::Point, 2> boundaryPoints;
	WGS84::Point meshOrigin;
	Eigen::Vector2i tileSize;
	bool useCoarseSrtmResolution;
	bool produceAsciiPly;
	int numThreads;
};

CommandLineArguments parseCommandLine(int argc, const char** argv);