
#ifdef _MSC_VER
#if _MSC_VER >= 1900
//On MSVC before 2015, docopt won't work
//(since the compiler doesn't support multiple features of C++11)
#define USE_DOCOPT
#endif
#endif

//Fallback solution for MSVC2013
//To be removed as soon as MSVC2015 turns out to be stable enough.
#ifndef USE_DOCOPT

#include "CLI.h"
#include "Parsing.h"

#include <vector>
#include <string>
#include <regex>
#include <iostream>
#include <bitset>

void printAndExit(const std::string &message)
{
	std::cout << message << std::endl;
	exit(0);
}

CommandLineArguments parseCommandLine(int argc, const char** argv)
{
	CommandLineArguments arguments;

	std::vector<std::string> args(argc, "");
	for (int i = 0; i < argc; ++i)
		args[i] = std::string(argv[i]);

	if (args.empty()) 
		printAndExit(USAGE);

	bool optionsParsed = false;
	//Agument[0] contains the executable path
	auto currentArgument = args.begin()+1;

	enum Option: unsigned int
	{
		HELP = 0,
		VERSON,
		ORIGIN,
		TILE_SIZE,
		ASCII_PLY,
		NUM_THREADS,
		COARSE_SRTM,
		
		FROM_POINT,
		TO_POINT,

		INPUT_DIRECTORY,
		OUTPUT_DIRECTORY,

		NUM_OPTIONS
	};

	std::bitset<NUM_OPTIONS> parsedOptions = 0;
	auto setOptionParsed =
		[&](Option option) {
			if (parsedOptions[option])
				printAndExit(USAGE);
			else
				parsedOptions[option] = true;
		};

	auto getArgumentValue =
		[&](const std::smatch &match)
		{
			if (match.str(1).empty())
				return *(++currentArgument);
			else
				return match.str(1);
		};

	for(currentArgument; currentArgument != args.end(); ++currentArgument)
	{
		static const auto &icase = std::regex_constants::icase;

		if (!optionsParsed)
		{
			std::smatch match;

			static const std::regex helpPattern("^-h|--help$", icase);
			if (std::regex_match(*currentArgument, helpPattern))
				printAndExit(USAGE);

			static const std::regex versionPattern("^-v|--version$", icase);
			if (std::regex_match(*currentArgument, versionPattern))
				printAndExit(VERSION);

			static const std::regex originPattern("^-o|(?:--origin=(\\S+))$");
			if (std::regex_match(*currentArgument, match, originPattern))
			{
				setOptionParsed(ORIGIN);
				auto argument = getArgumentValue(match);
				arguments.meshOrigin = parseWGS84(argument);
				continue;
			}

			static const std::regex tileSizePattern("^-s|(?:--tile-size=(\\S+))$");
			if (std::regex_match(*currentArgument, match, tileSizePattern))
			{
				setOptionParsed(TILE_SIZE);
				auto argument = getArgumentValue(match);
				arguments.tileSize = parseSize(argument);
				continue;
			}

			static const std::regex asciiPlyPattern("^-a|--ascii-ply$");
			if (std::regex_match(*currentArgument, match, asciiPlyPattern))
			{
				setOptionParsed(ASCII_PLY);
				arguments.produceAsciiPly = true;
				continue;
			}

			static const std::regex numThreadsPattern("^-j|(?:-j=(\\S+))$");
			if (std::regex_match(*currentArgument, match, numThreadsPattern))
			{
				setOptionParsed(NUM_THREADS);
				auto argument = getArgumentValue(match);
				arguments.numThreads = parseNumberOfThreads(argument);
				continue;
			}

			static const std::regex coarseSrtmPattern("^-srtm3$");
			if (std::regex_match(*currentArgument, match, coarseSrtmPattern))
			{
				setOptionParsed(COARSE_SRTM);
				arguments.useCoarseSrtmResolution = true;
				continue;
			}

			optionsParsed = true;
		}

		if (optionsParsed)
		{
			if (!parsedOptions[FROM_POINT])
				arguments.boundaryPoints[0] = parseWGS84(*currentArgument);
			else if (!parsedOptions[TO_POINT])
				arguments.boundaryPoints[1] = parseWGS84(*currentArgument);
			else if (!parsedOptions[INPUT_DIRECTORY])
				arguments.inputDirectory = *currentArgument;
			else if (!parsedOptions[OUTPUT_DIRECTORY])
				arguments.outputDirectory = *currentArgument;
			else //Too many options...
				printAndExit(USAGE);
		}
	}

	if (!parsedOptions[FROM_POINT] || !parsedOptions[TO_POINT])
		printAndExit(USAGE);

	auto &from = arguments.boundaryPoints[0];
	auto &to = arguments.boundaryPoints[1];
	if (!parsedOptions[ORIGIN])
		arguments.meshOrigin = WGS84::Point(.5*(from.latitude()  + to.latitude()),
		                                    .5*(from.longitude() + to.longitude()),
											.5*(from.altitude()  + to.altitude()));

	if (!parsedOptions[NUM_THREADS])
		arguments.numThreads = parseNumberOfThreads("#CPUs");

	return arguments;
}

#endif