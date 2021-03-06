#include <iostream>
#include <iomanip>
#include <NaiveMatrix.h>
#include <vector>
#include <algorithm>
#include <cmath>

using std::cout;
using std::cin;
using std::endl;
using std::boolalpha;
using std::vector;
using std::fill;
using std::copy;
using std::abs;

int main(int argc, char **argv)
{
    cout << std::setprecision(10);
    vector<long double> b;
    unsigned int n = 0, m = 0;
    cout << "n is number of rows, and m is number of columns" << endl;
    cout << "Enter n : ";
    cin >> n;
    cout << "Enger m : ";
    cin >> m;

    cout << "(n, m) is (" << n << "," << m << ")" << endl;
    NaiveMatrix<long double> nm(n, m);
    cout << "after : (n, m) is (" << n << "," << m << ")" << endl;
    b.clear();
    b.resize(n);

    for(unsigned int i = 0; i < n; ++i)
    {
        for(unsigned int j = 0; j < m; ++j)
        {
            cout << "enter element [" <<  i  << "," << j << "] = ";
            cin >> nm(i, j);
        }
    }
    nm.print();
    cout<< "norm is " << nm.norm1() << endl;

    for(int i = 0; i < n; ++i){
        cout << "enter element of b[" << i << "] = ";
        cin >> b[i];
    }

    nm.print();
    cout << "[";
    for(int i = 0; i < b.size(); ++i)
    {
        cout << b[i];
        if ((i+1) != b.size())
            cout << ","; 
    }
    cout << "]" << endl;

    bool checkDiagonallydominant = nm.checkDiagonallydominant();
    cout << "check to see if Matrix is diagonally dominant : " << boolalpha << checkDiagonallydominant << endl;

    if(!checkDiagonallydominant){
        return 0;
    }

    cout << "Which method do you want to solve the equation ? " << endl;
    cout << "1 - Jacobian method" << endl;    
    cout << "2 - Gauss–Seidel method" << endl;    
    int method = 0;
    cin >> method;

    if(method == 1)
    {
        // Jacobian method
        long double epsilon = 0;
        cout << "enter epsilon " << endl; 
        cin >> epsilon;
        cout << "value of epsilon set to " << epsilon << endl;
        vector<long double> answer;
        vector<long double> n_answer;
        answer.resize(n);
        fill(answer.begin(), answer.end(), 0);
        n_answer.resize(n);
        fill(n_answer.begin(), n_answer.end(), 0);
        bool end = false;
        while(end == false)
        {
            for(int i = 0; i < nm.get_row(); ++i)
            {
                n_answer[i] = b[i];
                for(int j = 0; j < nm.get_column(); ++j)
                {
                    if(i == j)
                        continue;
                    n_answer[i] += (- (nm(i, j) * answer[j]));
                }
                n_answer[i] /= nm(i, i);
            }
            long double max = 0;
            for(int i = 0; i < nm.get_row(); ++i)
            {
                if (max < abs(answer[i] - n_answer[i]))
                    max = abs(answer[i] - n_answer[i]);
            }
            if(epsilon > max)
                end = true;
            copy(n_answer.begin(), n_answer.end(), answer.begin());
            fill(n_answer.begin(), n_answer.end(), 0);
        }

        cout << "Jacobian method answers : ";
        for(int i = 0; i < answer.size(); ++i)
        {
            cout << answer[i];
            if ((i+1) != answer.size())
                cout << ",";
        }        
        cout << endl;
    }
    else if (method == 2)
    {
        // Gauss-Seidel method
        long double epsilon = 0;
        cout << "enter epsilon " << endl; 
        cin >> epsilon;
        cout << "value of epsilon set to " << epsilon << endl;
        vector<long double> answer, old_answer;
        answer.resize(n);
        old_answer.resize(n);
        fill(answer.begin(), answer.end(), 0);

        bool end = false;
        while(end == false)
        {
            copy(answer.begin(), answer.end(), old_answer.begin());
            for(int i = 0; i < nm.get_row(); ++i)
            {
                answer[i] = b[i];
                for(int j = 0; j < nm.get_column(); ++j)
                {
                    if(i == j)
                        continue;
                    answer[i] += (- (nm(i, j) * answer[j]));
                }
                answer[i] /= nm(i, i);
            }

            long double max = 0;
            for(int i = 0; i < nm.get_row(); ++i)
            {
                if (max < abs(answer[i] - old_answer[i]))
                    max = abs(answer[i] - old_answer[i]);
            }
            if(epsilon > max)
                end = true;
        }

        cout << "Gauss-Seidel method answers : ";
        for(int i = 0; i < answer.size(); ++i)
        {
            cout << answer[i];
            if ((i+1) != answer.size())
                cout << ",";
        }        
        cout << endl;
    }
    else{
        cout << "Invalid input, Fatal error" <<endl;
        return 0;
    }

    if(!nm.determinant()){
        cout << "determinant is zero, matrix is not inversable" << endl;
        exit(0);
    }
    
    NaiveMatrix<long double> nmi = nm;
    nmi.invert();
    cout << "inverse matrix is : " << endl;
    nmi.print();

    cout << "cond_1 is " << nm.norm1() * nmi.norm1() << endl;
    return 0;
}


