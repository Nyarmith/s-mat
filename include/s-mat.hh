#ifndef __S_MAT_HH__
#define __S_MAT_HH__

#include <type_traits>
#include <utility>

namespace SMat {
// our base abstraction is a the vector
template<typename T, unsigned Cols>
struct Vec
{
    constexpr Vec() {}

    explicit constexpr Vec(const T val)
    {
        for (unsigned i=0; i<Cols; ++i)
            data[i] = val;
    }

    constexpr Vec(const std::initializer_list<T> &l)
    {
        unsigned i=0;
        for (auto e : l)
        {
            if (i>= Cols)
                return;
            data[i] = e;
            ++i;
        }
    }

    constexpr T& operator[](unsigned i)
    {
        return data[i];
    }

    constexpr const T& operator[](unsigned i) const
    {
        return data[i];
    }

    constexpr T operator*(const Vec &o) const
    {
        T ret=0;
        for (unsigned i=0; i<Cols; ++i)
            ret += data[i]*o.data[i];
        return ret;
    }

    constexpr Vec operator/(const T o) const
    {
        Vec ret;
        for (unsigned i=0; i<Cols; ++i)
            ret[i] = data[i]/o;
        return ret;
    }

    constexpr Vec operator-() const
    {
        Vec ret;
        for (unsigned i=0; i<Cols; ++i)
            ret[i] = -data[i];
        return ret;
    }

    constexpr T size() const
    {
        T ret=0;
        for (unsigned i=0; i<ret; ++i)
            ret += data[i]*data[i];
        return sqrt(ret);
    }

    template<typename Fn,
    typename std::enable_if_t<std::is_invocable_r_v<T,Fn,T>>* = nullptr>
    constexpr Vec operator%(Fn f) const
    {
        Vec ret;
        for (unsigned i=0; i<Cols; ++i)
            ret[i] = f(data[i]);
        return ret;
    };

    T data[Cols]{0};
};

template<typename T, unsigned Cols>
constexpr Vec<T,Cols> operator*(const Vec<T,Cols> &A, const T s)
{
    Vec<T,Cols> ret;
    for (unsigned i=0; i<Cols; ++i)
        ret[i] = A.data[i]*s;
    return ret;
}

template<typename T, unsigned Cols>
constexpr Vec<T,Cols> operator*(const T s, const Vec<T,Cols> &A)
{
    return A*s;
}

template<typename T, unsigned Cols>
constexpr bool operator==(const Vec<T,Cols> &l, const Vec<T,Cols> &r)
{
    for (unsigned i=0; i<Cols; ++i)
        if (l[i] != r[i])
            return false;

    return true;
}

template<typename T, unsigned Cols>
constexpr bool operator!=(const Vec<T,Cols> &l, const Vec<T,Cols> &r)
{
    return !(l==r);
}

template<typename T, unsigned Cols>
constexpr Vec<T,Cols> operator+(const Vec<T,Cols> &l, const Vec<T,Cols> &r)
{
    Vec<T,Cols> ret;
    for (unsigned i=0; i<Cols; ++i)
        ret[i] = l[i] + r[i];
    return ret;
}

template<typename T, unsigned Cols>
constexpr Vec<T,Cols> operator-(const Vec<T,Cols> &l, const Vec<T,Cols> &r)
{
    Vec<T,Cols> ret;
    for (unsigned i=0; i<Cols; ++i)
        ret[i] = l[i] - r[i];
    return ret;
}

template<typename T, unsigned Rows, unsigned Cols=Rows>
struct Matrix
{
    Vec<T,Cols> mat[Rows];

    constexpr Matrix(){ }

    explicit constexpr Matrix (const T v)
    {
        for (int i=0; i<Rows; ++i)
        {
            mat[i] = Vec<T,Cols>(v);
        }
    }

    constexpr Matrix (const std::initializer_list<Vec<T,Cols>> &l)
    {
        unsigned r=0;
        for (const auto e : l)
        {
            if (r>=Rows)
                return;
            mat[r] = e;
            ++r;
        }
    }

    constexpr Matrix (const std::initializer_list<T> &l)
    {
        unsigned i=0;
        unsigned r=0;
        unsigned c=0;
        for (auto e : l)
        {
            if (i>=Cols*Rows)
                return;
            mat[r][c] = e;
            ++c;
            if (c >= Cols)
            {
               c=0;
               ++r;
            }
            ++i;
        }
    }

    constexpr Vec<T,Cols>& operator[](unsigned i)
    {
        return mat[i];
    }

    constexpr const Vec<T,Cols>& operator[](unsigned i) const
    {
        return mat[i];
    }

    constexpr const Vec<T,Cols>& row(int i) const
    {
        return mat[i];
    }

    constexpr Vec<T,Rows> col(int j) const
    {
        Vec<T,Rows> ret;
        for (unsigned r=0; r<Rows; ++r)
            ret[r] = mat[r][j];
        return ret;
    }

    constexpr Matrix<T,Cols,Rows> transpose() const
    {
        Matrix<T,Cols,Rows> ret;
        for (unsigned i=0; i<Cols; ++i)
            ret[i] = col(i);
    };

    template<typename Fn,
    typename std::enable_if_t<std::is_invocable_r_v<T,Fn,T>>* = nullptr>
    constexpr Matrix operator%(Fn f) const
    {
        Matrix ret;
        for (unsigned i=0; i<Rows; ++i)
            ret[i] = mat[i]%f;
        return ret;
    }
};

template <typename T, unsigned N, unsigned M, unsigned O>
constexpr Matrix<T,N,O> operator*(const Matrix<T,N,M> &A, const Matrix<T,M,O> &B)
{
    Matrix<T,N,O> ret;
    for (unsigned i=0; i<N; ++i)
        for (unsigned j=0; j<O; ++j)
            ret[i][j] = A[i]*B.col(j);
    return ret;
}

template <typename T, unsigned R, unsigned C>
constexpr Matrix<T,R,C> operator*(const Matrix<T,R,C> &A, const T s)
{
    Matrix<T,R,C> ret;
    for (unsigned i=0; i<R; ++i)
        ret[i] = A[i]*s;
    return ret;
}

template <typename T, unsigned R, unsigned C>
constexpr Matrix<T,R,C> operator*(const T s, const Matrix<T,R,C> &A)
{
    return A*s;
}

template <typename T, unsigned R, unsigned C>
constexpr Matrix<T,R,C> operator/(const Matrix<T,R,C> &A, const T s)
{
    Matrix<T,R,C> ret;
    for (unsigned i=0; i<R; ++i)
        ret[i] = A[i]/s;
    return ret;
}

template <typename T, unsigned R, unsigned C>
constexpr Matrix<T,R,C> operator+(const Matrix<T,R,C> &A, const Matrix<T,R,C> &B)
{
    Matrix<T,R,C> ret;
    for (unsigned i=0; i<C; ++i)
        ret[i] = A[i] + B[i];
    return ret;
}

template <typename T, unsigned R, unsigned C>
constexpr Matrix<T,R,C> operator-(const Matrix<T,R,C> &A, const Matrix<T,R,C> &B)
{
    Matrix<T,R,C> ret;
    for (unsigned i=0; i<C; ++i)
        ret[i] = A[i] - B[i];
    return ret;
}

template <typename T, unsigned R, unsigned C>
constexpr bool operator==(const Matrix<T,R,C> &A, const Matrix<T,R,C> &B)
{
    for (unsigned i=0; i<R; ++i)
        if (A[i] != B[i])
            return false;

    return true;
}

template <typename T, unsigned R, unsigned C>
constexpr bool operator!=(const Matrix<T,R,C> &A, const Matrix<T,R,C> &B)
{
    return !(A == B);
}

template <typename T, unsigned N, typename P,
typename std::enable_if_t<std::is_integral<P>::value>* = nullptr>
constexpr Matrix<T,N,N> pow(Matrix<T,N,N> A, P p)
{
    if (p == 0)
        return Matrix<T,N,N>(1);
    if (p == 1)
        return A;
    if (p == 2)
        return A*A;
    if (p%2 == 1)
        return A*(pow(A,p-1));
    else
        return pow(pow(A,p/2),2);
}

template <typename T, unsigned C>
constexpr bool is_close(Vec<T,C> A, Vec<T,C> B, const T delta)
{
    for (unsigned i=0; i<C; ++i)
        if (A[i] <= B[i] - delta || A[i] >= B[i] + delta)
            return false;
    return true;
}

template <typename T, unsigned R, unsigned C>
constexpr bool is_close(Matrix<T,R,C> A, Matrix<T,R,C> B, const T delta)
{
    for (unsigned i=0; i<R; ++i)
        if (!is_close(A[i],B[i],delta))
            return false;
    return true;
}

} // namespace SMat

#endif
