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
        Matrix(const std::vector<std::vector<double> >& values);

        // Destructor
        ~Matrix(){}

        // Operators
        std::vector<double>& operator[](size_t i);
        Matrix operator/(double val) const;
        Matrix operator+(const Matrix& rhs) const;
        Matrix operator+(double value) const;
        Matrix operator+(const std::vector<double> vect) const;
        Matrix operator-(const Matrix& rhs) const;
        Matrix operator-(double value) const;
        Matrix operator-() const;
        Matrix operator*(const Matrix& rhs) const;
        std::vector<double> operator*(const std::vector<double>& vect) const;
        Matrix operator*(double value) const;

        // Accessors
        size_t getNCols();
        size_t getNRows();

        // Mutators
        void setMatrix(const std::vector<std::vector<double> >& values);
        
        // Matrix Methods
        void display() const;
        Matrix invert() const;

    private:
        size_t m_nCols;
        size_t m_nRows;
        std::vector<std::vector<double> > m_M;
};