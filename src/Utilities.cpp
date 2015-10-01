
#include <exception>
#include <sstream>

#include "Utilities.h"

std::string absolutePath(const std::string &path)
{
    std::string result;
    char *abs = realpath(path.c_str(), NULL);
    if (!abs)
    {
        std::stringstream ss;
        ss << "Could not find file or directory \""
           << path << "\".";
        throw std::runtime_error(ss.str());
    }
    else
    {
        result = std::string(abs);
        free(abs);
    }
    return result;
}