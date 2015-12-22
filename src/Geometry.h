
#ifndef GEOMETRY_H
#define GEOMETRY_H

#if defined(_MSC_VER)
#define NOMINMAX
#endif

#include <algorithm>
#include <iterator>

#include <Eigen/Geometry>

/*! \file math/Geometry.h
 *
 *  \brief This file contains the necessary types to convert points between
 *         various geographic coordinate systems.
 *
 *  The respective coordinate systems are defined as classes that contain a
 *  specific nested vector implementation that is used to represent any point
 *  in their coordinate space.
 *
 *  As a central reference system, the Earth Centered, Earth Fixed cartesian
 *  coordinate system is used (cf. the EarthCenteredEarthFixed class).
 *  All coordinate system types have to implement the two methods toECEF(...)
 *  and fromECEF(...), which therefore allows the conversion between any two
 *  coordinate systems through the ECEF representation of the according point.
 *
 *  There are two basic types of coordinate systems:
 *
 *  - Statically defined coordinate systems such as ECEF and WGS84 usually follow
 *    a well-known definition and can therefore be accessed as Singletons through
 *    an according static instance() method. The conversion between points of two
 *    statically defined coordinate systems can be done implicitly.
 *
 *  - Dynamically defined coordinate systems depend on certain parameters and thus
 *    are only known at runtime. Thus, there are no possible implicit conversions
 *    of any point of such a coordinate system into any other representation.
 *    Conversions of points in those systems have therefore to be done explicitly
 *    through their respective ECEF representations.
 *
 *  Example code:
 *
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 *
 *  WGS84::Point sion(46.233333, 7.366654); //46.2333°N, 7.3667°E, 0m altitude
 *
 *  ECEF::Point p = sion; //Sion in it's ECEF representation
 *
 *  LCS localCoords(sion); //A local coordinate system with sion as it's origin
 *  LCS::Point p2 = localCoords.fromECEF(p);
 *
 *  //sion is approximately (46.2333 / 7.36665 / 0)
 *  //p    is approximately (4.38e+6 /  566678 / 4.58e+6)
 *  //p2   is approximately (0 / 0 / 0)
 *
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*! \brief A type that represents the geocentric cartesian coordinate system.
 *
 *  The ECEF coordinate system is used as the central reference system for all
 *  geometric points. Every other coordinate system has to implement the methods
 *  fromECEF(...) and toECEF(...) to allow a smooth conversion of points from one
 *  coordinate system to another.
 *
 *  The ECEF system is defined as the right-handed cartesian coordinate system
 *  in which the X axis passes through the intersection of the prime meridian and
 *  the equator, while the Z axis passes through the IERS reference pole (IRP).
 */
struct EarthCenteredEarthFixed
{
    /*! \brief A representation of a point in the ECEF space.
     *
     *  This Point class represents any point in the ECEF space through a
     *  position vector. The class provides an implicit conversion member
     *  for other point classes that are defined through static coordinate
     *  system classes.
     */
    struct Point: public Eigen::Vector3d
    {
        //! A type alias for the coordinate system of this type
        typedef EarthCenteredEarthFixed CoordinateSystem;

        Point() = default;
        Point(const Point &other) = default;
        Point(double x, double y, double z):
            Eigen::Vector3d(x,y,z)
        {}

        /*! \brief An explicit conversion method to create a point from any vector.
         *
         *  To avoid accidental "conversions" that only involve copying the coefficients of
         *  a vector, this method has to be called explicitly.
         */
        explicit Point(const Eigen::Vector3d &v):
            Eigen::Vector3d(v)
        {}

        Point &operator=(const Point &other) = default;

        //! Convert a point from another coordinate system to ECEF.
        template<class T>
        Point &operator=(const T &p)
        {
            (*this) = T::CoordinateSystem::instance().toECEF(p);
            return *this;
        }

        //! Convert a point from another coordinate system to ECEF.
        template<class T>
        Point(const T &p)
        {
            (*this) = p;
        }
    };

    //! The ECEF system is statically defined; one can access it's instance through a singleton function.
    static const EarthCenteredEarthFixed &instance();

    //! For ECEF, any conversions are simple identity transforms
    Point fromECEF(const EarthCenteredEarthFixed::Point &p) const { return p; }
    //! For ECEF, any conversions are simple identity transforms
    EarthCenteredEarthFixed::Point toECEF(const Point &p) const   { return p; }

private:
    EarthCenteredEarthFixed(){}
};
typedef EarthCenteredEarthFixed ECEF;

/*! \brief A type that represents the latest world geodetic system revision (WGS84, EPSG:4326).
 *
 *  WGS84 defines a standard spheroidal coordinate system, in which any point on earth can be
 *  described as a combination of latitude and longitude angles [°] and a altitude value [m] relative
 *  to a defined reference ellipsoid.
 */
struct WorldGeodeticSystem
{
    /*! \brief A representation of a point in the WGS84 coordinate system.
     *
     *  This Point class represents any point in the WGS84 space through a
     *  three dimensional vector containing the respective coordinate coefficients
     *  in the same order as defined in the ISO 6709 standard (latitude/longitude/altitude).
     */
    struct Point: public Eigen::Vector3d
    {
        //! A type alias for the coordinate system of this type
        typedef WorldGeodeticSystem CoordinateSystem;

        Point() = default;
        Point(const Point &other) = default;
        Point(double latitude, double longitude, double altitude = 0):
            Eigen::Vector3d(latitude,longitude,altitude)
        {}

        /*! \brief An explicit conversion method to create a point from any vector.
         *
         *  To avoid accidental "conversions" that only involve copying the coefficients of
         *  a vector, this method has to be called explicitly.
         */
        explicit Point(const Eigen::Vector3d &v):
            Eigen::Vector3d(v)
        {
            assert(v[0] >=  -90.0);
            assert(v[0] <=   90.0);
            assert(v[1] >= -180.0);
            assert(v[1] <=  180.0);
        }

        Point &operator=(const Point &other) = default;

        //! Convert a point from the another coordinate space to it's WGS84 representation.
        template<class T>
        Point &operator=(const T &p)
        {
            ECEF::Point ecef = p;
            (*this) = CoordinateSystem::instance().fromECEF(ecef);
            return *this;
        }

        //! Convert a point from the ECEF space to it's WGS84 representation.
        template<class T>
        Point(const T &p)
        {
            (*this) = p;
        }

        //! A simple access method for the first coefficient
        inline double &latitude()  { return (*this)[0]; }
        //! A simple access method for the second coefficient
        inline double &longitude() { return (*this)[1]; }
        //! A simple access method for the third coefficient
        inline double &altitude()    { return (*this)[2]; }

        //! A simple access method for the first coefficient
        inline double latitude()  const { return (*this)[0]; }
        //! A simple access method for the second coefficient
        inline double longitude() const { return (*this)[1]; }
        //! A simple access method for the third coefficient
        inline double altitude()    const { return (*this)[2]; }

        //! Check if the EGM-Point is valid (i.e. it's degree coefficients are in their defined range).
        bool isValid() const
        {
            //Check if latitude and longitude are well-formed (i.e. within their degree range)
            return std::abs(latitude()/90) <= 1.0 &&
                   std::abs(longitude()/180) <= 1.0;
        }
    };

    //! The WGS84 system is statically defined; one can access it's instance through a singleton function.
    static const WorldGeodeticSystem &instance();

    //! Convert the provided ECEF point into it's WGS84 representation.
    Point fromECEF(const ECEF::Point &ecef) const;
    //! Convert the provided WGS84 point into it's respective ECEF representation.
    ECEF::Point toECEF(const Point &p) const;

    /*! \brief Convert multiple ECEF points at once into their according WGS84 representations.
     *
     *  \arg srcFirst An iterator pointing to the first element to be transformed.
     *  \arg srcLast  An iterator pointing to the position after the last element to be transformed.
     *  \arg dstFirst An iterator pointing to the first element where the result should be written.
     */
    template<class I, class O>
    void fromECEF(I srcFirst, I srcLast, O dstFirst)
    {
        std::transform(srcFirst,srcLast,dstFirst,[&](decltype(*srcFirst) p){
            return fromECEF(ECEF::Point(p));
        });
    }

    /*! \brief Convert multiple ECEF points at once into their according WGS84 representations.
     *
     *  \arg src An iterable object containing all the points that should be transformed.
     *  \arg dst An iterable object of the same size where the according results are stored.
     */
    template<class I, class O>
    void fromECEF(const I &src, O &dst)
    {
        fromECEF(std::begin(src),std::end(src),std::begin(dst));
    }

    /*! \brief Convert multiple WGS84 points at once into their according ECEF representations.
     *
     *  \arg srcFirst An iterator pointing to the first element to be transformed.
     *  \arg srcLast  An iterator pointing to the position after the last element to be transformed.
     *  \arg dstFirst An iterator pointing to the first element where the result should be written.
     */
    template<class I, class O>
    void toECEF(I srcFirst, I srcLast, O dstFirst)
    {
        std::transform(srcFirst,srcLast,dstFirst,[&](decltype(*srcFirst) p){
            return toECEF(WorldGeodeticSystem::Point(p));
        });
    }

    /*! \brief Convert multiple WGS84 points at once into their according ECEF representations.
     *
     *  \arg src An iterable object containing all the points that should be transformed.
     *  \arg dst An iterable object of the same size where the according results are stored.
     */
    template<class I, class O>
    void toECEF(const I &src, O &dst)
    {
        toECEF(std::begin(src),std::end(src),std::begin(dst));
    }

private:
    WorldGeodeticSystem(){}
};
typedef WorldGeodeticSystem WGS84;

/*! \brief A type that represents the WGS84 coordinate system, but using altitude values based on the
*          EGM96 geoid.
*/
class GeoidConverter;
struct EarthGravitationalModel
{
    struct GeoidNotInitialized : public std::runtime_error
    {
        GeoidNotInitialized():
            std::runtime_error("")
        {}
    };

    /*! \brief A representation of a point in the WGS84-EGM96 coordinate system.
    *
    *  This Point class represents any point in the WGS84 space through a
    *  three dimensional vector containing the respective coordinate coefficients
    *  in the same order as defined in the ISO 6709 standard (latitude/longitude/altitude).
    */
    struct Point : public Eigen::Vector3d
    {
        //! A type alias for the coordinate system of this type
        typedef EarthGravitationalModel CoordinateSystem;

        Point() = default;
        Point(const Point &other) = default;
        Point(double latitude, double longitude, double altitude = 0) :
            Eigen::Vector3d(latitude, longitude, altitude)
        {}

        /*! \brief An explicit conversion method to create a point from any vector.
        *
        *  To avoid accidental "conversions" that only involve copying the coefficients of
        *  a vector, this method has to be called explicitly.
        */
        explicit Point(const Eigen::Vector3d &v):
            Eigen::Vector3d(v)
        {
            assert(v[0] >= -90.0);
            assert(v[0] <= 90.0);
            assert(v[1] >= -180.0);
            assert(v[1] <= 180.0);
        }

        Point &operator=(const Point &other) = default;

        //! Convert a point from the another coordinate space to it's EGM96 representation.
        template<class T>
        Point &operator=(const T &p)
        {
            ECEF::Point ecef = p;
            (*this) = CoordinateSystem::instance().fromECEF(ecef);
            return *this;
        }

        //! Convert a point from the another coordinate space to it's EGM96 representation.
        template<class T>
        Point(const T &p)
        {
            (*this) = p;
        }

        //! A simple access method for the first coefficient
        inline double &latitude()  { return (*this)[0]; }
        //! A simple access method for the second coefficient
        inline double &longitude() { return (*this)[1]; }
        //! A simple access method for the third coefficient
        inline double &altitude()    { return (*this)[2]; }

        //! A simple access method for the first coefficient
        inline double latitude()  const { return (*this)[0]; }
        //! A simple access method for the second coefficient
        inline double longitude() const { return (*this)[1]; }
        //! A simple access method for the third coefficient
        inline double altitude()    const { return (*this)[2]; }

        //! Check if the EGM-Point is valid (i.e. it's degree coefficients are in their defined range).
        bool isValid() const
        {
            //Check if latitude and longitude are well-formed (i.e. within their degree range)
            return std::abs(latitude()/90) <= 1.0 &&
                   std::abs(longitude()/180) <= 1.0;
        }
    };

    //! The WGS84 system is statically defined; one can access it's instance through a singleton function.
    static const EarthGravitationalModel &instance();

    //! Convert the provided ECEF point into it's EGM96 representation.
    Point fromECEF(const ECEF::Point &ecef) const;
    //! Convert the provided EGM96 point into it's respective ECEF representation.
    ECEF::Point toECEF(const Point &p) const;

    //! Convert the provided WGS84 point into it's EGM96 representation.
    Point fromWGS(const WGS84::Point &wgs) const;
    //! Convert the provided EGM96 point into it's respective WGS84 representation.
    WGS84::Point toWGS(const Point &p) const;

    /*! \brief Convert multiple ECEF points at once into their according WGS84 representations.
    *
    *  \arg srcFirst An iterator pointing to the first element to be transformed.
    *  \arg srcLast  An iterator pointing to the position after the last element to be transformed.
    *  \arg dstFirst An iterator pointing to the first element where the result should be written.
    */
    template<class I, class O>
    void fromECEF(I srcFirst, I srcLast, O dstFirst)
    {
        std::transform(srcFirst, srcLast, dstFirst, [&](decltype(*srcFirst) p){
            return fromECEF(ECEF::Point(p));
        });
    }

    /*! \brief Convert multiple ECEF points at once into their according WGS84 representations.
    *
    *  \arg src An iterable object containing all the points that should be transformed.
    *  \arg dst An iterable object of the same size where the according results are stored.
    */
    template<class I, class O>
    void fromECEF(const I &src, O &dst)
    {
        fromECEF(std::begin(src), std::end(src), std::begin(dst));
    }

    /*! \brief Convert multiple WGS84 points at once into their according ECEF representations.
    *
    *  \arg srcFirst An iterator pointing to the first element to be transformed.
    *  \arg srcLast  An iterator pointing to the position after the last element to be transformed.
    *  \arg dstFirst An iterator pointing to the first element where the result should be written.
    */
    template<class I, class O>
    void toECEF(I srcFirst, I srcLast, O dstFirst)
    {
        std::transform(srcFirst, srcLast, dstFirst, [&](decltype(*srcFirst) p){
            return toECEF(WorldGeodeticSystem::Point(p));
        });
    }

    /*! \brief Convert multiple WGS84 points at once into their according ECEF representations.
    *
    *  \arg src An iterable object containing all the points that should be transformed.
    *  \arg dst An iterable object of the same size where the according results are stored.
    */
    template<class I, class O>
    void toECEF(const I &src, O &dst)
    {
        toECEF(std::begin(src), std::end(src), std::begin(dst));
    }

        static bool initialized();
        static void loadGeoidData(const std::string &path = ".");
    private:
        EarthGravitationalModel();
        ~EarthGravitationalModel();

        static EarthGravitationalModel &mutableInstance();

        GeoidConverter *_geoidConverter;
};
typedef EarthGravitationalModel EGM96;

//! A template specialization that allows for a faster conversion between EGM96 and WGS84 points
template<>
inline EGM96::Point &EGM96::Point::operator= <WGS84::Point>(const WGS84::Point &wgs)
{
    (*this) = EGM96::instance().fromWGS(wgs);
    return *this;
}

//! A template specialization that allows for a faster conversion between EGM96 and WGS84 points
template<>
inline WGS84::Point &WGS84::Point::operator= <EGM96::Point>(const EGM96::Point &egm)
{
    (*this) = EGM96::instance().toWGS(egm);
    return *this;
}

/*! \brief A type that represents a local cartesian coordinate system that is defined by a single point on the geoid.
 *
 *  The local coordinate system is defined in such a way that it can be derived solely from a specific
 *  position on the geoid: The specified point acts as the new origin, with the Z axis pointing in the
 *  opposite direction of the center of mass of the earth, while the X and Y axes point in the east
 *  and north direction on the resulting plane.
 *
 *  The mapping of coordinates from the ECEF to the local space thus can be done through a simple
 *  affine 3D transformation (consisting of a translation and two rotations).
 */
struct LocalCoordinateSystem
{
    /*! \brief A representation of a point in a local coordinate system.
     *
     *  Note that this Point class only represents a point for a specific LCS instance
     *  through a position vector. The according LCS instance has to be known through
     *  the context of the point itself, as there are no stored references in this type.
     *
     *  Since the conversion of instances of this class may vary on runtime, there are
     *  no methods for implicit conversions between the according LCS and the ECEF space.
     *  Those conversions have to be done explicitly through the fromECEF(...) and toECEF(...)
     *  methods of the coordinate system instance itself.
     */
    struct Point: public Eigen::Vector3d
    {
        //! A type alias for the coordinate system of this type
        typedef LocalCoordinateSystem CoordinateSystem;

        Point() = default;
        Point(const Point &other) = default;
        Point(double x, double y, double z):
            Eigen::Vector3d(x,y,z)
        {}

        /*! \brief An explicit conversion method to create a point from any vector.
         *
         *  To avoid accidental "conversions" that only involve copying the coefficients of
         *  a vector, this method has to be called explicitly.
         */
        explicit Point(const Eigen::Vector3d &v):
            Eigen::Vector3d(v)
        {}

        Point &operator=(const Point &other) = default;
    };

    //! Construct a default LCS instance (located at 0°N, 0°E, 0m altitude)
    LocalCoordinateSystem():
        LocalCoordinateSystem(WGS84::Point(0,0,0))
    {}
    //! Construct a LCS instance corresponding to the provided WGS84 position.
    LocalCoordinateSystem(const WGS84::Point &origin);
    LocalCoordinateSystem(const LocalCoordinateSystem &other) = default;
    LocalCoordinateSystem &operator=(const LocalCoordinateSystem &other) = default;

    inline const WGS84::Point &origin() const { return _origin; }

    //! Convert the provided ECEF point into it's according LCS representation (for this instance).
    Point fromECEF(const ECEF::Point &p) const;
    //! Convert the provided LCS point into it's respective ECEF representation.
    ECEF::Point toECEF(const Point &p) const;

    /*! \brief Convert multiple ECEF points at once into their according LCS representations (for this instance).
     *
     *  \arg srcFirst An iterator pointing to the first element to be transformed.
     *  \arg srcLast  An iterator pointing to the position after the last element to be transformed.
     *  \arg dstFirst An iterator pointing to the first element where the result should be written.
     */
    template<class I, class O>
    void fromECEF(I srcFirst, I srcLast, O dstFirst)
    {
        std::transform(srcFirst,srcLast,dstFirst,[&](const ECEF::Point &p){
            return Point(_transform*p);
        });
    }

    /*! \brief Convert multiple ECEF points at once into their according LCS representations (for this instance).
     *
     *  \arg src An iterable object containing all the points that should be transformed.
     *  \arg dst An iterable object of the same size where the according results are stored.
     */
    template<class I, class O>
    void fromECEF(const I &src, O &dst)
    {
        fromECEF(std::begin(src),std::end(src),std::begin(dst));
    }

    /*! \brief Convert multiple LCS points at once into their according ECEF representations.
     *
     *  \arg srcFirst An iterator pointing to the first element to be transformed.
     *  \arg srcLast  An iterator pointing to the position after the last element to be transformed.
     *  \arg dstFirst An iterator pointing to the first element where the result should be written.
     */
    template<class I, class O>
    void toECEF(I srcFirst, I srcLast, O dstFirst)
    {
        auto transform = _transform.inverse();
        std::transform(srcFirst,srcLast,dstFirst,[&](const Point &p){
            return ECEF::Point(transform*p);
        });
    }

    /*! \brief Convert multiple LCS points at once into their according ECEF representations.
     *
     *  \arg src An iterable object containing all the points that should be transformed.
     *  \arg dst An iterable object of the same size where the according results are stored.
     */
    template<class I, class O>
    void toECEF(const I &src, O &dst)
    {
        toECEF(std::begin(src),std::end(src),std::begin(dst));
    }

    /*! \brief Convert multiple points from the current LCS to another one.
     *
     *  \arg srcFirst An iterator pointing to the first element to be transformed.
     *  \arg srcLast  An iterator pointing to the position after the last element to be transformed.
     *  \arg dstFirst An iterator pointing to the first element where the result should be written.
     *
     *  Note that this method is significantly faster than the equivalent conversions that result from
     *  the conversion to ECEF and to the target LCS from there, since the affine 3D transforms for
     *  both systems can be combined into one step.
     */
    template<class I, class O>
    void convertTo(const LocalCoordinateSystem &system, I srcFirst, I srcLast, O dstFirst)
    {
        Eigen::Affine3d transform = _transform.inverse()*system._transform;
        std::transform(srcFirst,srcLast,dstFirst,[&](decltype(*srcFirst) p){
            return Point(transform*Point(p));
        });
    }

    /*! \brief Convert multiple points from the current LCS to another one.
     *
     *  \arg src An iterable object containing all the points that should be transformed.
     *  \arg dst An iterable object of the same size where the according results are stored.
     *
     *  Note that this method is significantly faster than the equivalent conversions that result from
     *  the conversion to ECEF and to the target LCS from there, since the affine 3D transforms for
     *  both systems can be combined into one step.
     */
    template<class I, class O>
    void convertTo(const LocalCoordinateSystem &system, const I &src, O &dst)
    {
        convertTo(system,std::begin(src),std::end(src),std::begin(dst));
    }

private:
    WGS84::Point _origin;
    Eigen::Affine3d _transform;
};
typedef LocalCoordinateSystem LCS;

#endif //GEOMETRY_H
