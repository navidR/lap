#include <vector>
#include <numeric>
#include <functional>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iterator>


using std::vector;
using std::copy;
using std::cout;
using std::cin;
using std::endl;
using std::abs;
using std::exit;
using std::multiplies;

template <typename Element>
class NaiveMatrix{
    private:
        // n: row
        // m: column
        const int n,m;
        std::vector<Element> internal;
        unsigned int index(unsigned int, unsigned int) const;
        bool isDiagonallydominant;
        const Element twoXtwoDeterminant(Element, Element, Element, Element) const;

    public:
        NaiveMatrix(const NaiveMatrix&);
        NaiveMatrix(const int, const int);
        Element& operator()(const int, const int);
        void print();
        const int get_row() const;
        const int get_column() const;
        const bool checkDiagonallydominant() const;
        const Element determinant() const;
        void invert();
        void multiplyToScalar(Element);
        void minor();
        void cofactor();
        void adjugate();
        void transpose();
        Element norm1();
};

#include <NaiveMatrix.cpp>