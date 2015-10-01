#include "Geometry.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include <GeographicLib/GeoCoords.hpp>
#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/Geoid.hpp>

const WorldGeodeticSystem &WorldGeodeticSystem::instance()
{
    static WorldGeodeticSystem instance;
    return instance;
}

const EarthCenteredEarthFixed &EarthCenteredEarthFixed::instance()
{
    static EarthCenteredEarthFixed instance;
    return instance;
}

const EarthGravitationalModel &EarthGravitationalModel::instance()
{
    return mutableInstance();
}

EarthGravitationalModel &EarthGravitationalModel::mutableInstance()
{
    static EarthGravitationalModel instance;
    return instance;
}

WGS84::Point WorldGeodeticSystem::fromECEF(const ECEF::Point &ecef) const
{
    Point wgs;
    GeographicLib::Geocentric::WGS84().Reverse(ecef.x(), ecef.y(), ecef.z(),
                                               wgs.latitude(), wgs.longitude(), wgs.altitude());
    return wgs;
}

ECEF::Point WorldGeodeticSystem::toECEF(const Point &wgs) const
{
    EarthCenteredEarthFixed::Point ecef;
    GeographicLib::Geocentric::WGS84().Forward(wgs.latitude(), wgs.longitude(), wgs.altitude(),
                                               ecef.x(), ecef.y(), ecef.z());
    return ecef;
}


class GeoidConverter
{
public:
    GeoidConverter(const std::string &path = ".") :
        geoid("egm96-5", path, true, true)
    {}

    GeographicLib::Geoid geoid;
};


EarthGravitationalModel::EarthGravitationalModel() :
_geoidConverter(NULL)
{}

EarthGravitationalModel::~EarthGravitationalModel()
{
    if (_geoidConverter)
        delete _geoidConverter;
}

bool EarthGravitationalModel::initialized()
{
    return instance()._geoidConverter;
}

void EarthGravitationalModel::loadGeoidData(const std::string &path)
{
    auto instance = EarthGravitationalModel::mutableInstance();

    if (!instance._geoidConverter)
        instance._geoidConverter = new GeoidConverter(path);

    instance._geoidConverter->geoid.CacheAll();
}

EGM96::Point EarthGravitationalModel::fromECEF(const ECEF::Point &ecef) const
{
    return fromWGS(ecef);
}

ECEF::Point EarthGravitationalModel::toECEF(const Point &egm) const
{
    return toWGS(egm);
}


EarthGravitationalModel::Point EarthGravitationalModel::fromWGS(const WGS84::Point &wgs) const
{
    assert(initialized());

    auto egm = Point(wgs);
    egm.altitude() = _geoidConverter->geoid.ConvertHeight(wgs.latitude(), wgs.longitude(), wgs.altitude(), GeographicLib::Geoid::ELLIPSOIDTOGEOID);
    return egm;
}

WGS84::Point EarthGravitationalModel::toWGS(const EarthGravitationalModel::Point &p) const
{
    assert(initialized());

    auto wgs = WGS84::Point(p);
    wgs.altitude() = _geoidConverter->geoid.ConvertHeight(wgs.latitude(), wgs.longitude(), wgs.altitude(), GeographicLib::Geoid::GEOIDTOELLIPSOID);
    return wgs;
}


LocalCoordinateSystem::LocalCoordinateSystem(const WGS84::Point &origin):
    _origin(origin)
{
    WGS84::Point north(origin); north.latitude() += 1.0;

    ECEF::Point center (origin);
    WGS84::Point up = origin; up.altitude() += 100;
    Eigen::Vector3d normal = (ECEF::Point(up)-center).normalized();
    Eigen::Vector3d y = (ECEF::Point(north)-center).normalized(); //Approximative north direction (w/ wrong inclination)
    Eigen::Vector3d x = y.cross(normal).normalized();             //East direction

    //Correct the north vector to be aligned with the plane defined by our origin and normal
    y = Eigen::AngleAxisd(-(M_PI/2-acos(y.dot(normal))),x)*y;

    _transform = Eigen::Translation3d(-center);

    //Rotate the mesh to align the northward tangent with the Y axis in the local coordinate system
    Eigen::Vector3d transformed = _transform*(center+y);
    _transform = Eigen::AngleAxisd(acos(Eigen::Vector3d::UnitY().dot(transformed)), transformed.cross(Eigen::Vector3d::UnitY()).normalized()) * _transform;

    //Then, align the eastward tangent with the X axis
    transformed = (_transform*(center+x)).normalized();
    _transform = Eigen::AngleAxisd(acos(Eigen::Vector3d::UnitX().dot(transformed)), transformed.cross(Eigen::Vector3d::UnitX()).normalized()) * _transform;
}

LCS::Point LocalCoordinateSystem::fromECEF(const ECEF::Point &p) const
{
    return Point(_transform*p);
}

ECEF::Point LocalCoordinateSystem::toECEF(const Point &p) const
{
    return ECEF::Point(_transform.inverse()*p);
}
