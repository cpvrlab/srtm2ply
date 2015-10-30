
#ifndef MESH_H
#define MESH_H

#if defined(_MSC_VER)
#define NOMINMAX
#endif

#include "Geometry.h"
#include "SRTM.h"

#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <array>

#include "Utilities.h"

struct Mesh
{
    typedef Eigen::Vector3f Point;
    typedef std::array<size_t,3> Face;

    typedef std::vector<Point> VertexPositions;
    typedef std::vector<Face>  Faces;

    VertexPositions vertexPositions;
    Faces faces;

    Mesh() = default;
    Mesh(const Mesh&) = default;
	Mesh(Mesh &&other):
		vertexPositions(std::move(other.vertexPositions)),
		faces(std::move(other.faces))
	{}

    Mesh &operator=(const Mesh&) = default;
	Mesh &operator=(Mesh &&other)
	{
		vertexPositions = std::move(other.vertexPositions);
		faces = std::move(other.faces);
	}

    void toAsciiPly(std::ostream &output) const;
    void toBinaryPly(std::ostream &output) const;

private:
    //Determine the float precision for the ASCII ply format
    int determinePrecision() const;
    size_t determineIndexSize() const;
};

#endif //MESH_H
