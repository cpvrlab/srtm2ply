
#ifdef _MSC_VER
#if _MSC_VER >= 1900
//On MSVC before 2015, docopt won't work
//(since the compiler doesn't support multiple features of C++11)
#define USE_DOCOPT
#endif
#endif

#ifdef USE_DOCOPT

#include "CLI.h"
#include "Parsing.h"

#include <docopt.h>

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

			for (auto &str : list)
				std::cout << "                        " << str << ",\n";
			std::cout << "                  ]\n\n";
		}
		else if (!value) std::cout << "       (empty)\n\n";
		else std::cout << "       (unknown)\n\n";

	}
	std::cout << std::endl;
}

CommandLineArguments parseCommandLine(int argc, const char** argv)
{
	CommandLineArguments arguments;

	auto args = docopt::docopt(USAGE, { argv + 1, argv + argc },
								true,                       // show help if requested
								VERSION                     // version string
							  );

	auto cwd = getCwd();

	arguments.inputDirectory  = args["<INPUT_DIRECTORY>" ].isString() ? args["<INPUT_DIRECTORY>" ].asString() : cwd;
	arguments.outputDirectory = args["<OUTPUT_DIRECTORY>"].isString() ? args["<OUTPUT_DIRECTORY>"].asString() : cwd;

	arguments.inputDirectory  = absolutePath(arguments.inputDirectory);
	arguments.outputDirectory = absolutePath(arguments.outputDirectory);

	WGS84::Point from, to, origin;
	arguments.boundaryPoints[0] = parseWGS84(args["<p1>"].asString());
	arguments.boundaryPoints[1] = parseWGS84(args["<p2>"].asString());

	if (args["<origin>"].isString())
		arguments.meshOrigin = parseWGS84(args["<origin>"].asString());
	else
		arguments.meshOrigin = WGS84::Point(Eigen::Vector3d((from + to) / 2.0));

	if (args["--tile-size"].isString())
		arguments.tileSize = parseSize(args["--tile-size"].asString());
	else
		arguments.tileSize = Eigen::Vector2i(-1,-1);

	arguments.tileSize = parseNumberOfThreads(args["-j"].asString());
	
	arguments.useCoarseSrtmResolution = args["--srtm3"] && args["--srtm3"].asBool();
	arguments.produceAsciiPly = args["--ascii-ply"] && args["--ascii-ply"].asBool();
}

#endif
