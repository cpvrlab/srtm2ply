
#include <climits>
#include <cstdlib>
#include <exception>
#include <sstream>

#include "Utilities.h"

#if defined(_MSC_VER)
#include <windows.h>
#include <Shlwapi.h>
#include <vector>
#endif

std::string getAbsolutePath(const std::string &path)
{
	std::string result = "";
#if defined(_MSC_VER)
	result.resize(MAX_PATH,' ');
	SetLastError(0);
	GetFullPathName(path.c_str(), result.size(), &result[0], nullptr);
	if (GetLastError())
		throw std::runtime_error("Could not determine absolute path of file.");
	if (!PathFileExists(result.c_str()))
		result = "";
#else
	char *abs = realpath(path.c_str(), NULL);
	if (abs)
	{
		result = std::string(abs);
		free(abs);
	}
#endif
	return result;
}

std::string absolutePath(const std::string &path)
{
    std::string result = getAbsolutePath(path);
	if (result.empty())
	{
		std::stringstream ss;
		ss << "Could not find file or directory \""
			<< path << "\".";
		throw std::runtime_error(ss.str());
	}
	return result;
}