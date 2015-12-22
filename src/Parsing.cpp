
#include "Parsing.h"

#include <regex>
#include <sstream>
#include <cmath>
#include <ctgmath>
#include <thread>

std::string parseErrorMessage(const std::string &expression)
{
    std::stringstream ss;
    ss << "Cannot parse expression \"" << expression << "\".";
    return ss.str();
}

ParseError::ParseError(const std::string &expression):
    std::runtime_error(parseErrorMessage(expression))
{}

WGS84::Point parseWGS84(const std::string &coordinates)
{
    std::regex pattern(
                    //Latitude
                    //--------
                       "^"
                       "=?"
                       "("
                           "-?"                    //Sign
                           "(?:"
                               "\\d+(?:\\.\\d+)?"  //Decimal (i.e. 1, 19.2, 0.002...)
                               "|"
                               "\\.\\d+"           //Short-hand small decimal (.002...)
                           ")"
                       ")"
                       "([NS])"                    //Latitude direction (north, south)

                       "\\s*"                      //Possible whitespace

                     //Longitude
                     //---------
                       "("
                           "-?"                    //Sign
                           "(?:"
                               "\\d+(?:\\.\\d+)?"  //Decimal (i.e. 1, 19.2, 0.002...)
                               "|"
                               "\\.\\d+"           //Short-hand small decimal (.002...)
                           ")"
                       ")"
                       "([EW])"                    //Longitude direction (east, west)

                       "\\s*"                      //Possible whitespace

                     //Altitude (optional)
                     //-------------------
                       "(?:"
                           "("
                               "-?"                    //Sign
                               "(?:"
                                   "\\d+(?:\\.\\d+)?"  //Decimal (i.e. 1, 19.2, 0.002...)
                                   "|"
                                   "\\.\\d+"           //Short-hand small decimal (.002...)
                               ")"
                           ")"
                           "m"                         //Mandatory unit specification (Meters)
                       ")?"
                       "$"
                       );

    std::smatch match;
    if (!std::regex_match(coordinates, match, pattern))
        throw ParseError(coordinates);

    double latitude = strtod(std::string(match[1]).c_str(),nullptr);
    char latitudeDirection = std::string(match[2])[0];

    double longitude = strtod(std::string(match[3]).c_str(),nullptr);
    char longitudeDirection = std::string(match[4])[0];

    double altitude = strtod(std::string(match[5]).c_str(),nullptr);

    if (!std::isfinite(latitude) || abs(latitude) > 180)
        throw ParseError(coordinates);

    if (!std::isfinite(longitude) || abs(longitude) > 90)
        throw ParseError(coordinates);

    if (!std::isfinite(altitude))
        throw ParseError(coordinates);

    latitude  *= latitudeDirection  == 'N' ? 1 : -1;
    longitude *= longitudeDirection == 'E' ? 1 : -1;

    return WGS84::Point(latitude,longitude,altitude);
}


Eigen::Vector2i parseSize(const std::string &size)
{
    std::regex pattern( "^"
                        "=?"
                        "("
                            "\\d+"
                        ")"
                        "x"
                        "("
                            "\\d+"
                        ")"
                        "$"
                     );

    std::smatch match;
    if (!std::regex_match(size, match, pattern))
        throw ParseError(size);

    int x,y;
    x = strtol(std::string(match[1]).c_str(),nullptr,10);
    y = strtol(std::string(match[2]).c_str(),nullptr,10);

    if (x <= 0 || y <= 0)
        throw ParseError(size);

    return {x,y};
}

int parseNumberOfThreads(const std::string &n)
{
    std::regex pattern1("^=?#CPUs$");
    std::regex pattern2("^=?(\\d+)$");

    std::smatch match;

    if (std::regex_match(n, match, pattern1))
        return std::thread::hardware_concurrency();

    if (!std::regex_match(n, match, pattern2))
        return 1;

    int result = strtol(std::string(match[1]).c_str(),nullptr,10);
    return std::max(result,1);
}

