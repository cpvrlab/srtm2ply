
#include "Mesh.h"


int Mesh::determinePrecision() const
{
    Eigen::AlignedBox3f bounds;
    for (auto &pos: vertexPositions)
        bounds.extend(pos);

    if (bounds.isEmpty()) return 0;
    float minWidth = bounds.diagonal().minCoeff();
    int precision = 0;
    float threshold = 1000.0f;
    while(minWidth < threshold)
    {
        precision++;
        threshold /= 10;
    }
    return precision;
}

size_t Mesh::determineIndexSize() const
{
    size_t numVertices = vertexPositions.size();
    for (size_t i = 1; i < sizeof(size_t); i *= 2)
    {
        size_t threshold = 1;
        threshold <<= 8*i;
        if (numVertices <= threshold)
            return i;
    }
    return sizeof(size_t);
}


//Create a binary faces buffer with indexes of the provided type (e.g. uint8_t)
template<class T>
std::vector<uint8_t> createBinaryFacesBufferForType(const Mesh::Faces &faces)
{
    const auto indexSize = sizeof(T);
    const auto faceSize = sizeof(Mesh::Face)/sizeof(typename Mesh::Face::value_type);
    const auto bufferElementSize = faceSize*indexSize+1;

    std::vector<uint8_t> buffer(faces.size()*bufferElementSize);
    int i = 0;
    for (auto &face: faces)
    {
        uint8_t *ptr = &buffer[i++ * bufferElementSize];

        uint8_t &size = ptr[0];
        size = (uint8_t) face.size();

        for (unsigned int j = 0; j < face.size(); ++j)
        {
            T &index = *((T*)(ptr+1+j*indexSize));
            index = (T) face[j];
        }
    }
    return buffer;
}

template<class... T>
inline std::vector<unsigned char> createBinaryFacesBuffer(size_t indexSize, const Mesh::Faces &faces);

//Iterate through the template arguments, check the size of the type
//and generate a buffer if it matches
template<class T, class... U>
inline std::vector<unsigned char> _createBinaryFacesBuffer(size_t indexSize, const Mesh::Faces &faces)
{
    if (indexSize == sizeof(T))
        return createBinaryFacesBufferForType<T>(faces);
    else
        return createBinaryFacesBuffer<U...>(indexSize,faces);
}

template<class... T>
inline std::vector<unsigned char> createBinaryFacesBuffer(size_t indexSize, const Mesh::Faces &faces)
{
    return _createBinaryFacesBuffer<T...>(indexSize,faces);
}

//There was no matching index type provided, throw an exception
template<>
inline std::vector<unsigned char> createBinaryFacesBuffer<>(size_t, const Mesh::Faces&)
{
    throw std::runtime_error("Unknown index size.");
}


void Mesh::toAsciiPly(std::ostream &output) const
{
    output << "ply\n";
    output << "format ascii 1.0\n";
    output << "element vertex " << vertexPositions.size() << '\n';
    output << "property float32 x\n";
    output << "property float32 y\n";
    output << "property float32 z\n";
    output << "element face " << faces.size() << '\n';
    output << "property list uint8 int32 vertex_indices\n";
    output << "end_header\n";

    int precision = determinePrecision();
    for (auto &v: vertexPositions)
    {
        for (int i = 0; i < 3; ++i)
            output << std::setprecision(precision) << std::showpoint << std::fixed << v[i] << ' ';
        output << '\n';
    }

    for (auto &face: faces)
    {
        output << face.size() << ' ';
        for (int i = 0; i < 3; ++i)
            output << face[i] << ' ';
        output << '\n';
    }
    output << '\n';
}

void Mesh::toBinaryPly(std::ostream &output) const
{
    auto indexSize = determineIndexSize();

    output << "ply\n";
    output << "format "
           << (Architecture::ENDIANNESS == Architecture::E_BIG_ENDIAN ? "binary_big_endian" : "binary_little_endian")
           << " 1.0\n";
    output << "element vertex " << vertexPositions.size() << '\n';
    output << "property float32 x\n";
    output << "property float32 y\n";
    output << "property float32 z\n";
    output << "element face " << faces.size() << '\n';
    output << "property list uint8 uint" << (8*indexSize) << " vertex_indices\n";
    output << "end_header\n";

    output.write((const char*)vertexPositions.data(),vertexPositions.size()*sizeof(decltype(vertexPositions.front())));

    std::vector<unsigned char> facesBuffer = createBinaryFacesBuffer<uint8_t, uint16_t, uint32_t, uint64_t>(indexSize,faces);
    output.write((const char*)facesBuffer.data(), facesBuffer.size()*sizeof(decltype(facesBuffer.front())));
}
