#include <vector>
#include <iostream>
#include <cmath>


using std::vector;
using std::cout;
using std::cin;
using std::endl;
using std::abs;

template <typename Element>
class NaiveMatrix{
    private:
        const int n,m;
        std::vector<Element> internal;
        unsigned int index(unsigned int, unsigned int) const;
        bool isDiagonallydominant;
    public:
        NaiveMatrix(const int, const int);
        Element& operator()(const int, const int);
        void print();
        const int get_row() const;
        const int get_column() const;
        const bool checkDiagonallydominant() const;
};

#include <NaiveMatrix.cpp>