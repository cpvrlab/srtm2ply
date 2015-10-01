
#ifndef ARCHITECTURE_H
#define ARCHITECTURE_H

#include <cstdint>

namespace Architecture
{
    enum Endianness: uint8_t
    {
        E_LITTLE_ENDIAN = 0,
        E_BIG_ENDIAN
    };

    inline Endianness determineEndianness()
    {
        short i = 1;
        return (Endianness)(*((unsigned char*)&i) == 0);
    }

    static const Endianness ENDIANNESS = determineEndianness();
}

#endif //SRTM_H
