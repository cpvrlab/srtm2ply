
#ifndef PARSING_H
#define PARSING_H

#if defined(_MSC_VER)
#define NOMINMAX
#endif

#include <exception>

#include "Geometry.h"

struct ParseError: public std::runtime_error
{
    ParseError(const std::string &coordinates);
};

WGS84::Point parseWGS84(const std::string &coordinates);

Eigen::Vector2i parseSize(const std::string &size);

int parseNumberOfThreads(const std::string &n);

#endif // PARSING_H
