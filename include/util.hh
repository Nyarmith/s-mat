#ifndef __S_MAT_UTIL_HH__
#define __S_MAT_UTIL_HH__

#include "s-mat.hh"
#include <iostream>


namespace smat {

template<typename T, unsigned C>
std::ostream& operator<<(std::ostream& os, const Vec<T,C> o)
{
    for (unsigned c=0; c<C; ++c)
    {
        if (c != 0) os << ", ";
        os << o[c] << " ";
    }
    os << "\n";
    return os;
}

template<typename T, unsigned R, unsigned C>
std::ostream& operator<<(std::ostream& os, const Matrix<T,R,C> o)
{
    for (unsigned r=0; r<R; ++r)
        os << o[r];
    return os;
}
} // namespace smat

#endif
