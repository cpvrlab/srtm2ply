
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

	bool parsing = true;
	bool optionsParsed = false;
	std::string *arg = &args.front();
	while (parsing)
	{
		static const auto &icase = std::regex_constants::icase;

		if (!optionsParsed)
		{
			std::smatch match;

			auto &currentArgument = *arg;

			static const std::regex helpPattern("^-h|--help$", icase);
			if (std::regex_match(currentArgument, helpPattern))
				printAndExit(USAGE);

			static const std::regex versionPattern("^-v|--version$", icase);
			if (std::regex_match(currentArgument, versionPattern))
				printAndExit(VERSION);

			static const std::regex originPattern("^-o|(?:--origin=.*)$");
			if (std::regex_match(currentArgument, match, originPattern))
			{
				
			}
		}

		if (optionsParsed)
		{

		}
	}

	return arguments;
}

#endif