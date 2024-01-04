#pragma once
#include <cstddef>
#include <vector>

class Matrix
{
    public:
        // Default Constructor
        Matrix();

        // Initializer Constructors
        Matrix(size_t nRows, size_t nCols);
        Matrix(std::vector<std::vector<double> >& values);

        // Destructor
        virtual ~Matrix(){}

        // Operators
        std::vector<double>& operator[](size_t i);
        Matrix operator/(const double val);
        Matrix operator+(Matrix& rhs);
        Matrix operator+(const double value);
        Matrix operator+(std::vector<double> vect);
        Matrix operator-(Matrix& rhs);
        Matrix operator-(const double value);
        Matrix operator-();
        Matrix operator*(Matrix& rhs);
        std::vector<double> operator*(std::vector<double>& vect);
        Matrix operator*(const double value);

        // Accessors
        double getNCols();
        double getNRows();

        // Mutators
        void setMatrix(std::vector<std::vector<double> >& values);
        
        // Matrix Methods
        void display();
        double determinant();
        Matrix cofactorMatrix();
        Matrix transpose();
        Matrix invert();

    private:
        size_t m_nCols;
        size_t m_nRows;
        std::vector<std::vector<double> > m_M;
};