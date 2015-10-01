
#ifndef SRTM_H
#define SRTM_H

#include "Geometry.h"
#include "Architecture.h"

#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <memory>
#include <mutex>
#include <iostream>
#include <iomanip>
#include <thread>
#include <mutex>
#include <memory>
#include <exception>

#include "Tensor.h"

template<class T, Architecture::Endianness E, int R>
class SrtmTile
{
public:
    typedef Eigen::Vector2i SamplingPointLocation;
    typedef Eigen::AlignedBox2i Bounds;
    typedef T Element;
    typedef SrtmTile<T,E,R> Tile;
    typedef Tensor<T,2> Data;
    typedef typename Data::Size Size;
    typedef typename Data::Index Index;
    typedef typename Data::Position Position;

    static const Architecture::Endianness SOURCE_ENDIANNES = E;
    static const size_t RESOLUTION = R;

    //Exceptions
    //----------

    //! A runtime error signalling problems with the file system
    struct FileNotFound: public std::runtime_error
    {
        FileNotFound(const std::string &file):
            std::runtime_error(message(file))
        {}

        static std::string message(const std::string &file)
        {
            std::stringstream ss;
            ss << "Could not load SRTM terrain file \""
               << file << "\": The specified file doesn't exist.";
            return ss.str();
        }
    };

    //! A runtime error signalling a corrupt data file
    struct InvalidData: public std::runtime_error
    {
        InvalidData(const std::string &what):
            std::runtime_error(what)
        {}
    };

    template<class B, class U>
    class IteratorBase : public std::iterator < std::bidirectional_iterator_tag, T, std::ptrdiff_t, T*, T& >
    {
    public:
        typedef B *Block;
        typedef U (B::*IndexAccessMethod)(Index) const;
        typedef U Element;

        IteratorBase() :
            _block(nullptr),
            _index(0),
            _accessMethod(nullptr)
        {}

        IteratorBase(Block block, Index index, IndexAccessMethod method) :
            _block(block),
            _index(index),
            _accessMethod(method)
        {}

        IteratorBase(const IteratorBase&) = default;
        IteratorBase &operator=(const IteratorBase&) = default;

        bool operator==(const IteratorBase &other) const
        {
            return _block == other._block
                && _index == other._index
                && _accessMethod == other._accessMethod;
        }
        bool operator!=(const IteratorBase &other) const { return !(*this == other); }

        Element operator->() { return (_block->*_accessMethod)(_index); }
        Element operator*()  { return (_block->*_accessMethod)(_index); }

        IteratorBase &operator ++ ()
        {
            ++_index;
            return *this;
        }
        IteratorBase operator ++ (int) { IteratorBase tmp(*this); ++(*this); return tmp; }
        IteratorBase &operator -- ()
        {
            --_index;
            return *this;
        }
        IteratorBase operator -- (int) { IteratorBase tmp(*this); ++(*this); return tmp; }

        void swap(IteratorBase &other) { IteratorBase tmp(*this); *this = other; other = tmp; }

    private:
        Block _block;
        Index _index;
        IndexAccessMethod _accessMethod;
    };

    template<class U>
    class RangeBase
    {
    public:
        typedef U (Tile::*IndexAccessMethod)(Index) const;
        typedef IteratorBase<const Tile, U> Iterator;

        RangeBase(Tile &tile, IndexAccessMethod iAccess):
            _tile(&tile),
            _indexAccess(iAccess)
        {}

        Iterator begin() const { return Iterator(_tile, 0,                  _indexAccess); }
        Iterator end()   const { return Iterator(_tile, _tile->numValues(), _indexAccess); }

        Iterator cbegin() const { return begin(); }
        Iterator cend()   const { return end();   }

        const T &operator[](Index index)         const { assert(_tile); return (_tile->*_indexAccess)(index);  }
        const T &operator[](const Position &pos) const { assert(_tile); return (_tile->*_indexAccess)(_tile->indexAtPosition(index)); }

    private:
        Tile *_tile;
        IndexAccessMethod    _indexAccess;
    };

    struct AltitudeRange: public RangeBase<const Element&>
    {
        typedef RangeBase<const Element&> Base;

        AltitudeRange(Tile &tile):
            Base(tile, &Tile::altitudeAt)
        {}
    };

    struct WgsRange: public RangeBase<WGS84::Point>
    {
        typedef RangeBase<WGS84::Point> Base;

        WgsRange(Tile &tile):
            Base(tile, &Tile::wgsAt<true>)
        {}
    };

    class Subview
    {
    public:
        typedef Tile::Bounds Bounds;
        typedef Eigen::Vector2i Size;
        typedef Eigen::Vector2i Offset;

        Subview(const Tile &tile, const Bounds &bounds):
            _tile(&tile),
            _bounds(_tile->_bounds.intersection(bounds)),
            _size( _bounds.isEmpty() ? Size(0,0)
                                     : _bounds.max()-_bounds.min()+Eigen::Vector2i(1,1)),
            _offset(_bounds.min() - _tile->_bounds.min())
        {}

        Subview(const Tile &tile, const WGS84::Point &from, const WGS84::Point &to):
            Subview(Tile::bounds(from,to))
        {}

        Subview(const Subview&) = default;
        Subview &operator=(const Subview&) = default;

        inline bool isEmpty() const noexcept { return _bounds.isEmpty(); }
        inline const Bounds &bounds() const noexcept { return _bounds; }
        inline const Offset &offset() const noexcept { return _offset; }
        inline const Size &size() const noexcept { return _size; }
        inline size_t numValues() const noexcept { return _size.prod(); }
        inline Index indexAtPosition(const Position &pos) const noexcept { assertCorrectPosition(pos); return Data::positionToIndex(pos,_size); }
        inline Position positionAtIndex(Index index) const noexcept { assertCorrectIndex(index); return Data::indexToPosition(index,_size); }
        inline bool contains(const WGS84::Point &p) const { return Tile::Contains(_bounds,p); }

        Position tilePosition(const Position &pos) const
        {
            assertCorrectPosition(pos);
            return pos + Eigen::Array2i(_offset);
        }
        inline Index tileIndex(const Position &pos) const { return _tile->indexAtPosition(tilePosition(pos)); }
        inline Index tileIndex(Index index) const { return tileIndex(positionAtIndex(index)); }
        inline Position tilePosition(Index index) const { return tilePosition(positionAtIndex(index)); }

        template<bool FETCH_HEIGHT = true>
        inline WGS84::Point wgsAt(Index index)  const { assertCorrectIndex(index); return _tile->wgsAt<FETCH_HEIGHT>(tilePosition(index)); }

        template<bool FETCH_HEIGHT = true>
        WGS84::Point wgsAt(const Position &pos) const { assertCorrectPosition(pos); return _tile->wgsAt<FETCH_HEIGHT>(tilePosition(pos)); }

        template<bool FETCH_HEIGHT = true>
        inline WGS84::Point southwest() const { return wgsAt<FETCH_HEIGHT>(Position{0,0}); }

        template<bool FETCH_HEIGHT = true>
        inline WGS84::Point southeast() const { return wgsAt<FETCH_HEIGHT>(Position{0, size()[1]-1}); }

        template<bool FETCH_HEIGHT = true>
        inline WGS84::Point northwest() const { return wgsAt<FETCH_HEIGHT>(Position{size()[0]-1, 0}); }

        template<bool FETCH_HEIGHT = true>
        inline WGS84::Point northeast() const { return wgsAt<FETCH_HEIGHT>(Position{size()[0]-1, size()[1]-1}); }

        inline const T &altitudeAt(Index index) const { return altitudeAt(positionAtIndex(index)); }
        inline const T &altitudeAt(const Position &pos) const { return _tile->altitudeAt(tilePosition(pos)); }
        inline double altitudeAt(const WGS84::Point &p) const { assert(contains(p)); return _tile->altitudeAt(p); }
        inline void setToGroundLevel(WGS84::Point &p) const { p.altitude() = altitudeAt(p); }

    private:
        inline void assertCorrectPosition(const Position &position) const
        {
            assert((position >= 0).all());
            assert((position < Eigen::Array2i(_size)).all());
        }

        inline void assertCorrectIndex(Index index) const
        {
            assert(index >= 0);
            assert(index < numValues());
        }

        const Tile *_tile;
        Bounds _bounds;
        Size _size;
        Offset _offset;
    };

    class Cache
    {
    public:
        typedef typename std::shared_ptr<Tile> TilePtr;
        typedef typename std::weak_ptr<Tile> WeakTilePtr;

        Cache(std::string path):
            _path(std::move(path))
        {}

        TilePtr getTile(const SamplingPointLocation &p)
        {
            Eigen::Vector2i pos = Tile::fileOrigin(p);

            std::lock_guard<std::mutex> lock(_mutex);
            auto &ptr = _cache[pos];
            auto tile = ptr.lock();
            if (!tile)
            {
                if (_blacklist[pos])
                    return nullptr;

                try
                {
                    tile = std::make_shared<Tile>(Tile::loadFromDisk(p,_path));
                    ptr = tile;
                } catch(FileNotFound &e)
                {
                    _blacklist[pos] = true;
                    return nullptr;
                }
            }
            return tile;
        }

    private:
        //! Implement a hashing function for tile identifiers
        struct Hasher {
           size_t operator()(const Eigen::Array2i &p) const {
              std::hash<int> hash;
              auto h1 = hash(p[0]);
              auto h2 = hash(p[1]);
              //Imitate hash_combine() from boost to combine hashes from both coefficients
              return h1 + 0x9e3779b9 + (h2 << 6) + (h2 >> 2);
           }

           size_t operator()(const Eigen::Vector2i &p) const {
              return this->operator()(Eigen::Array2i(p));
           }
        };

        std::mutex _mutex;
        std::string _path;
        std::unordered_map<Eigen::Vector2i,WeakTilePtr,Hasher> _cache;
        std::unordered_map<Eigen::Vector2i,bool,Hasher> _blacklist;
    };

    //Static methods
    //--------------

    static SrtmTile loadFromDisk(const WGS84::Point &position, const std::string &path = ".")
    {
        return loadFromDisk(fileOrigin(position),path);
    }

    static SrtmTile loadFromDisk(const Eigen::Vector2i &origin, const std::string &path = ".")
    {
        std::string file_location = fileName(origin,path);
        std::ifstream file(file_location,std::ifstream::binary);

        //Check that the file exists
        if (!file)
            throw FileNotFound(file_location);

        //Check if the file has the right format
        file.seekg(0, std::ios::end);
        size_t fileSize = file.tellg();
        if (fileSize != sizeof(Element)*R*R)
            throw InvalidData("Error while parsing SRTM terrain file: The file size doesn't match the specified format.");

        //Allocate the tile in memory
        SrtmTile tile(origin, {R, R});

        //Jump to the start of the file and read it into memory

        //This would be faster but read in the data flipped (along the Y axis):
        file.seekg(0, std::ios::beg);
        file.read((char*)&tile._data[0],tile._data.numValues()*sizeof(Element));

        //Correct the endianness if necessary
        tile.correctEndianness();

        //While the defined origin for each SRTM tile is generally the south west corner
        //(i.e. the lower left corner), the "pixel origin" of the height data contained
        //in the files is in the top left corner (i.e. north west).
        //Thus, we have to flip the image along the Y axis to allow simpler position
        //calculations further on.
        tile.flipYAxis();

        return tile;
    }

    static SrtmTile stitch(const Bounds &bounds, Cache &cache)
    {
        SrtmTile result(bounds);
        auto tiles = loadIntersectingTiles(bounds, cache);

        for (auto &tile: tiles)
            if (tile)
                copyIntersection(*tile,result);

        return result;
    }

    static SrtmTile stitch(const Bounds &bounds, const std::string &path)
    {
        Cache cache(path);
        return stitch(bounds,cache);
    }

    static SrtmTile stitch(const WGS84::Point &from, const WGS84::Point &to, Cache &cache)
    {
        return stitch(bounds(from,to), cache);
    }

    static SrtmTile stitch(const WGS84::Point &from, const WGS84::Point &to, const std::string &path)
    {
        Cache cache(path);
        return stitch(from,to,cache);
    }

    static SamplingPointLocation nearestSamplingPoint(const WGS84::Point &p)
    {
        static const auto round = [](double a) { return std::round(a); };
        return samplingPoint(p,round);
    }

    static SamplingPointLocation southwestSamplingPoint(const WGS84::Point &p)
    {
        static const auto floor = [](double a) { return std::floor(a); };
        return samplingPoint(p,floor);
    }

    static SamplingPointLocation northeastSamplingPoint(const WGS84::Point &p)
    {
        static const auto ceil  = [](double a) { return std::ceil(a);  };
        return samplingPoint(p,ceil);
    }

    static WGS84::Point coordinate(const Eigen::Vector2i &samplingPoint)
    {
        Eigen::Vector2d coords = samplingPoint.cast<double>() / (R-1);
        return WGS84::Point(coords[1],coords[0],0);
    }

    static Bounds bounds(const SamplingPointLocation &southwest, const Size &size)
    {
        return Bounds(southwest).extend(southwest+Eigen::Vector2i(size)-Eigen::Vector2i(1,1));
    }

    static Bounds bounds(const WGS84::Point &from, const WGS84::Point &to)
    {
        WGS84::Point max, min;
        min.head(2) = from.head(2).cwiseMin(to.head(2));
        max.head(2) = from.head(2).cwiseMax(to.head(2));

        auto sw = southwestSamplingPoint(min);
        auto ne = northeastSamplingPoint(max);
        return bounds(sw, ne-sw+Eigen::Vector2i(1,1));
    }

    SrtmTile(const SrtmTile&) = default;
    SrtmTile(SrtmTile&&) = default;

    SrtmTile &operator=(const SrtmTile&) = default;
    SrtmTile &operator=(SrtmTile&&) = default;

    inline bool isEmpty() const noexcept { return _bounds.isEmpty(); }
    inline const Bounds &bounds() const noexcept { return _bounds; }
    inline const Size &size() const noexcept { return _data.size(); }
    inline size_t numValues() const noexcept { return _data.numValues(); }
    inline Index indexAtPosition(const Position &pos) const noexcept { return _data.indexAtPosition(pos); }
    inline Position positionAtIndex(Index index) const noexcept { return _data.positionAtIndex(index); }
    inline bool contains(const WGS84::Point &p) const { return contains(_bounds,p); }

    template<bool FETCH_HEIGHT = true>
    inline WGS84::Point wgsAt(Index index) const { return wgsAt<FETCH_HEIGHT>(positionAtIndex(index)); }

    template<bool FETCH_HEIGHT = true>
    WGS84::Point wgsAt(const Position &pos) const
    {
        assert((pos >= 0).all());
        assert((pos < size()).all());

        auto wgs = coordinate(_bounds.min()+Eigen::Vector2i(pos));
        if (FETCH_HEIGHT)
            wgs.altitude() = altitudeAt(pos);
        return wgs;
    }

    template<bool FETCH_HEIGHT = true>
    inline WGS84::Point southwest() const { return wgsAt<FETCH_HEIGHT>(Position{0,0}); }

    template<bool FETCH_HEIGHT = true>
    inline WGS84::Point southeast() const { return wgsAt<FETCH_HEIGHT>(Position{0, size()[1]-1}); }

    template<bool FETCH_HEIGHT = true>
    inline WGS84::Point northwest() const { return wgsAt<FETCH_HEIGHT>(Position{size()[0]-1, 0}); }

    template<bool FETCH_HEIGHT = true>
    inline WGS84::Point northeast() const { return wgsAt<FETCH_HEIGHT>(Position{size()[0]-1, size()[1]-1}); }

    const T &altitudeAt(Index index) const
    {
        assert(index >= 0);
        assert(index < numValues());
        return _data[index];
    }

    const T &altitudeAt(const Position &pos) const
    {
        assert((pos >= 0).all());
        assert((pos < size()).all());
        return _data[pos];
    }

    /*! \brief Obtain the altitude value of a specific WGS84 location.
     *
     *  The altitude value is calculated through a bilinear interpolation
     *  over it's neighbor samples.
     */
    inline double altitudeAt(const WGS84::Point &p) const
    {
        static const auto floor = [](double a) { return std::floor(a); };

        Eigen::Vector2d coords = (R-1)*p.head(2);
        std::swap(coords[0],coords[1]); //Invert latitude and longitude
        Eigen::Vector2d sw = coords.unaryExpr(floor);
        Eigen::Array2d offsets = coords-sw;
        Eigen::Array2i pos = sw.cast<int>() - _bounds.min();

        assert((pos >= 0).all());
        assert((size > pos).all());

        if (!(offsets > 0.0).any())
            return altitudeAt(pos);

        assert((size-pos > 1).all());

        return    offsets[0] *   offsets[1] *altitudeAt(Position(pos[0]  ,pos[1]  ))
             + (1-offsets[0])*   offsets[1] *altitudeAt(Position(pos[0]+1,pos[1]  ))
             +    offsets[0] *(1-offsets[1])*altitudeAt(Position(pos[0]  ,pos[1]+1))
             + (1-offsets[0])*(1-offsets[1])*altitudeAt(Position(pos[0]+1,pos[1]+1));
    }

    //! Correct (overwrite) the altitude value of the provided WGS84 point.
    inline void setToGroundLevel(WGS84::Point &p) const { p.altitude() = altitudeAt(p); }

    const AltitudeRange &altitudes() const { return _aRange; }
    const WgsRange &coordinates() const { return _wRange; }

private:
    static bool contains(const Eigen::AlignedBox2i &bounds, const WGS84::Point &p)
    {
        auto sw = southwestSamplingPoint(p);
        auto ne = northeastSamplingPoint(p);
        return bounds.contains(sw) && bounds.contains(ne);
    }

    static void copyIntersection(const SrtmTile &src, SrtmTile &dst)
    {
        auto bounds = src._bounds.intersection(dst._bounds);

        if (bounds.isEmpty())
            return;

        Tile::Subview srcView(src, bounds);
        Tile::Subview dstView(dst, bounds);

        const int n = srcView.numValues();
        for (int i = 0; i < n; ++i)
            dst._data[dstView.tileIndex(i)] = srcView.altitudeAt(i);
    }

    static Tensor<std::shared_ptr<Tile>,2> loadIntersectingTiles(const Bounds &bounds, Cache &cache)
    {
        Tensor<std::shared_ptr<Tile>,2> tiles;
        Eigen::Vector2i diagonal = bounds.max() - fileOrigin(bounds.min());
        Size size = diagonal / R + Eigen::Vector2i(1,1);
        tiles.assign(size,nullptr);

        for (int x = 0; x < size[0]; ++x)
        {
            for (int y = 0; y < size[1]; ++y)
            {
                Eigen::Vector2i origin = bounds.min()+Eigen::Vector2i((R-1)*x,(R-1)*y);
                origin = fileOrigin(origin);
                tiles[{x,y}] = cache.getTile(origin);
            }
        }
        return tiles;
    }

    //! Generate the SRTM file name for a specific WGS84 point.
    static std::string fileName(const Eigen::Vector2i &p, const std::string &path)
    {
        Eigen::Vector2i origin = fileOrigin(p) / (R-1);
        std::swap(origin[0],origin[1]); //Invert latitude and longitude
        std::stringstream ss;

        ss << path << "/";

        ss << (origin[0] >= 0? 'N' : 'S') << std::setw(2) << std::setfill('0') << abs(origin[0]);
        ss << (origin[1] >= 0? 'E' : 'W') << std::setw(3) << std::setfill('0') << abs(origin[1]);

        ss << ".hgt";
        return ss.str();
    }

    template<class F>
    static SamplingPointLocation samplingPoint(const WGS84::Point &p, F round)
    {
        Eigen::Vector2d coords(p[1],p[0]); //Invert latitude and longitude
        return Eigen::Vector2d((coords * (R-1)).unaryExpr(round)).cast<int>();
    }

    static Eigen::Vector2i fileOrigin(const WGS84::Point &p)
    {
        return fileOrigin(nearestSamplingPoint(p));
    }

    static Eigen::Vector2i fileOrigin(const Eigen::Vector2i &p)
    {
        auto ifloor = [](int a) {
            int offset = a % (R-1);
            return a - (offset + (bool(a < 0)*(R-1)));
        };
        return p.unaryExpr(ifloor);
    }

    //! Determine if the read data neets swapping.
    inline static bool isSwapNecessary() noexcept
    {
        auto endianness = Architecture::ENDIANNESS;
        return endianness != SOURCE_ENDIANNES;
    }

    void flipYAxis()
    {
        auto size = _data.size();
        std::vector<T> buffer(size[0],0);
        for (int y = 0; y < size[1]/2; ++y)
        {
            //Swapping elements manually (slow as hell):
//            for (int x = 0; x < size[0]; ++x)
//                std::swap(_data[{x,y}],_data[{x,size[1]-(y+1)}]);

            //Swap the data row-wise (mush faster)
            memcpy(&buffer[0]             , &_data[{0,y}]          , sizeof(T)*size[0]);
            memcpy(&_data[{0,y}]          , &_data[{0,size[1]-y-1}], sizeof(T)*size[0]);
            memcpy(&_data[{0,size[1]-y-1}], &buffer[0]             , sizeof(T)*size[0]);
        }
    }

    //! Switch the byte order of the sample values if necessary
    void correctEndianness()
    {
        static const bool SWAP_BYTES = isSwapNecessary();
        if (SWAP_BYTES)
        {
            std::transform(_data.cbegin(), _data.cend(), _data.begin(), [](const T &value) {
                T tmp;
                for (int i = 0; i < int(sizeof(T)/2); ++i)
                {
                    int j = sizeof(T)-1-i;
                    *((char*)&tmp + i) = *((char*)&value + j);
                    *((char*)&tmp + j) = *((char*)&value + i);
                }
                return tmp;
            });
        }
    }

    SrtmTile(const Bounds &bounds):
        _bounds(bounds),
        _data(),
        _aRange(*this),
        _wRange(*this)
    {
        Size size = bounds.diagonal()+Eigen::Vector2i(1,1);
        assert((size > 0).all());
        _data.assign(size,0);
    }

    SrtmTile(const Eigen::Vector2i &origin, const Size &size):
        SrtmTile(Bounds(origin,origin+Eigen::Vector2i(size[0]-1,size[1]-1)))
    {}

    Eigen::AlignedBox2i _bounds;
    Tensor<T,2> _data;

    AltitudeRange _aRange;
    WgsRange _wRange;
};

namespace SRTM1
{
    typedef SrtmTile<short,Architecture::E_BIG_ENDIAN,3601> Tile;
}

namespace SRTM3
{
    typedef SrtmTile<short,Architecture::E_BIG_ENDIAN,1201> Tile;
}

#endif //SRTM_H
