

/*
n is number of row
m is number of column
[0 1 2 3 
 4 5 6 7
 8 9 a b
 c d e f]
 
 turns to :
 vector of : [0, 1, 2, 4, 5, 6, 7, 8, 9, a, b, c, d, e, f];

*/

template <typename Element>
NaiveMatrix<Element>::NaiveMatrix(const int _n, const int _m) : n(_n), m(_m), isDiagonallydominant(false)
{
    internal.resize(n * m);
};

template<typename Element>
NaiveMatrix<Element>::NaiveMatrix(const NaiveMatrix& _n) : n( _n.get_row()), m(_n.get_column())
{
    internal = _n.internal;
};

template <typename Element>
Element& NaiveMatrix<Element>::operator()(const int in, const int im)
{
    if(in >= n || im >= m)
    {
        throw 0;
    }
    return internal[index(in, im)];
}

template <typename Element>
void NaiveMatrix<Element>::print()
{
    cout << "[";
    for(int i = 0; i < n; ++i)
    {
        if (i != 0)
            cout << " " ;
        for(int j = 0; j < m; ++j)
        {
            cout << internal[index(i, j)] ;
            if ((j + 1) != m)
                cout << ", ";
        }
        if ((i + 1) != n)
            cout << endl;
    }
    cout << "]" << endl;
}

template <typename Element>
unsigned int NaiveMatrix<Element>::index(unsigned int i, unsigned int j) const
{
    return (i * this->m) + j;
}

template <typename Element>
const int NaiveMatrix<Element>::get_row() const
{
    return this->n;
}

template <typename Element>
const int NaiveMatrix<Element>::get_column() const
{
    return this->m;
}

template <typename Element>
const Element NaiveMatrix<Element>::twoXtwoDeterminant(Element a, Element b, Element c, Element d) const
{
        return (a*d) - (b*c);
}


template <typename Element>
const bool NaiveMatrix<Element>::checkDiagonallydominant() const
{
    for(int i = 0; i < n; ++i)
    {
        Element sum = 0;
        for(int j = 0; j < m; ++j)
        {
            if(j == i)
                continue;
            sum += abs(internal[index(i, j)]);
        }
        if(abs(internal[index(i, i)]) <= sum)
            return false;
    }
    return true;
}

template <typename Element>
const Element NaiveMatrix<Element>::determinant() const
{
    if(this->m != this->n)
    {
        cout << "determinant only defined for square matrices, fatal error" << endl; 
        exit(0);       
    }
    else if(this->n > 3)
    {
        cout << "I am not capable of computing determinant of matrix with dimention bigger than 3, fatal error" << endl; 
        exit(0);       
    }
    Element det = 0;

    if(this->n == 2)
    {
        Element a = internal[index(0, 0)];
        Element b = internal[index(0, 1)];
        Element c = internal[index(1, 0)];
        Element d = internal[index(1, 1)];
        twoXtwoDeterminant(a, b, c, d);
    }

    // I did use sarrus-way of finding determinant of matrix with dimention 3
    vector<Element> sarrus;
    sarrus.resize(this->m);

    for(int step = 0; step < this->m; ++step)
    {
        fill(sarrus.begin(), sarrus.end(), 0);
        for(int i = 0, s = 0; i < this->n; ++i, ++s)
        {
            sarrus[s] = internal[index(i, (i + step) % this->m )];
        }
        det += accumulate(sarrus.begin(), sarrus.end(), 1, multiplies<Element>());
    }

    for(int step = 0; step < this->m; ++step)
    {
        fill(sarrus.begin(), sarrus.end(), 0);
        for(int i = (this->n - 1), j = 0, s = 0; i >= 0 && j < this->m; --i, ++j, ++s)
        {
          sarrus[s] = internal[index(i, (j + step) % this->m)];
        }
        det -= accumulate(sarrus.begin(), sarrus.end(), 1, multiplies<Element>());
    }
    return det;
}

template <typename Element>
void NaiveMatrix<Element>::minor()
{
    vector<Element> minors;
    minors.resize(internal.size());
    vector<Element> adbc;
    adbc.resize(4);
    for(int i = 0; i < this->n; ++i)
    {
        for(int j = 0; j < this->m; ++j)
        {
            int _index = 0;
            fill(adbc.begin(), adbc.end(), 0);
            for(int ii = 0; ii < this->n; ++ii)
            {   
                if(i == ii)
                    continue;
                for(int jj = 0; jj < this->m; ++jj)
                {
                    if(j == jj)
                        continue;
                    adbc[_index] = internal[index(ii, jj)];
                    _index++;
                }
            }
            minors[index(i, j)] = twoXtwoDeterminant(adbc[0], adbc[1], adbc[2], adbc[3]);
        }
    }
    internal = minors;
}

template <typename Element>
void NaiveMatrix<Element>::cofactor()
{
    for(int i = 0; i < internal.size(); ++i)
        internal[i] = i % 2 == 0 ? internal[i] * (+1) : internal[i] * (-1);
}

template <typename Element>
void NaiveMatrix<Element>::adjugate()
{
    this->transpose();
}

template <typename Element>
void NaiveMatrix<Element>::transpose()
{
    vector<Element> transposed;
    transposed.resize(internal.size());
    for(int i = 0; i < this->n; ++i)
    {
        for(int j = 0; j < this->m; ++j)
        {
            transposed[index(j, i)] = internal[index(i ,j)];
        }
    }
    internal = transposed;
}

template <typename Element>
void NaiveMatrix<Element>::invert()
{
    // copy to new matrix
    // minor
    // cofactors
    // adjugate
    // multiplyToScalar( determinant )
    if(this->m != this->n)
    {
        cout << "fatal error, incapable of finding invert of non-square matrix" << endl;
        exit(0);
    }
    Element det = this->determinant();
    if(!det)
    {
        cout << "determinant is Zero, there is no invert for this Matrix" << endl;
        exit(0);
    }
    
    if(this->m == 2)
    {
        this->multiplyToScalar(1 /(Element) det);
        return;
    }

    // 3x3 matrix 
    this->minor();
    this->cofactor();
    this->adjugate();
    this->multiplyToScalar(1 / det);
    return;
}

template <typename Element>
void NaiveMatrix<Element>::multiplyToScalar(Element e)
{
    // I think vs15 is not C++11, I am no sure 100%
    //std::for_each(internal.begin(), internal.end(), [](int &n){ n *= e;});
    //for(auto &n : internal)
    //    n *= e;
    for(int i = 0; i < internal.size(); ++i)
        internal[i] *= e;
}

template <typename Element>
Element NaiveMatrix<Element>::norm1()
{
    Element max = 0;
    for(int j = 0; j < this->m; ++j)
    {
        Element t_max = 0;
        for(int i = 0; i < this->n; ++i)
        {
            t_max += abs(internal[index(i, j)]);
        }
        if(t_max > max)
            max = t_max;
    }
    return max;
}