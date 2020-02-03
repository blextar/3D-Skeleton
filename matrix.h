#ifndef MATRIX_H
#define MATRIX_H

#include<iostream>
#include<vector>
#include<iterator>
#include<algorithm>

template<typename T> class matrix;
namespace std { template<typename T> void swap(matrix<T>&, matrix<T>&); };

template<typename T>
class matrix
{
    std::vector<T> buf;
    unsigned height_, width_;
    typedef T& ref;
    typedef T const & const_ref;
    public:
    matrix(): buf(), height_(0), width_(0) {}
    matrix(unsigned h, unsigned w, T val=T()): buf(h*w,val), height_(h), width_(w) {}

    void resize(unsigned h, unsigned w, T val=T())
    {   //TODO reorder data so that resize is not destructive
        buf.resize(h*w,val);
        height_=h;
        width_=w;
    };

    unsigned height() const {return height_; }
    unsigned width() const {return width_; }
    T* operator [](unsigned i) { return &(buf[i*width_]); }
    T const* operator [](unsigned i) const { return &(buf[i*width_]); }

    T& operator () (unsigned i, unsigned j) { if (j>=width_ || i>=height_) ((char*)0)[0]=1; return buf[i*width_+j]; }
    T const& operator () (unsigned i, unsigned j) const { if (j>=width_ || i>=height_) ((char*)0)[0]=1; return buf[i*width_+j]; }

    friend void std::swap<T>(matrix<T>&, matrix<T>&);
};

template<typename T>
std::ostream& operator << (std::ostream& out, matrix<T> const& M)
{
    const unsigned h=M.height(), w=M.width();
    out << h << ' ' << w << '\n';
    for (unsigned i=0; i!=h; ++i)
    {
        std::copy(M[i],M[i]+w,std::ostream_iterator<T>(out," "));
        out << "\n";
    }
    return out;
}

namespace std
{
    template<typename T>
    void swap(matrix<T>& a, matrix<T>& b)
    {
        swap(a.buf, b.buf);
        swap(a.height_ ,b.height_);
        swap(a.width_ ,b.width_);
    }
};

#endif // MATRIX_H
