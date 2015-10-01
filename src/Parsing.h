
#ifndef PARSING_H
#define PARSING_H

#include <exception>

#include "Geometry.h"

struct ParseError: public std::runtime_error
{
    ParseError(const std::string &coordinates);
};

WGS84::Point parseWGS84(const std::string &coordinates);

Eigen::Vector2i parseSize(const std::string &size);

#endif // PARSING_H
