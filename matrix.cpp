#include <iostream>
#include <vector>
#include "Matrix.h"

using namespace std;

// Default constructor
Matrix::Matrix() : m_nCols(3), m_nRows(3), m_M({{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}) { }

// Initializer Constructors
Matrix::Matrix(size_t nRows, size_t nCols) : m_nCols(nCols), m_nRows(nRows)
{
    m_M.resize(m_nRows, vector<double>(m_nCols, 0.0));
}

Matrix::Matrix(vector<vector<double> >& values) : m_nCols(values[0].size()), m_nRows(values.size()), m_M(values) { }

// Operators
vector<double>& Matrix::operator[](size_t i) { return this->m_M[i]; }

Matrix Matrix::operator/(const double val)
{
    Matrix res(this->m_nRows, this->m_nCols);

    for(size_t i = 0; i < this->m_nRows; i++)
    {
        for(size_t j = 0; j < this->m_nCols; j++)
        {
            res.m_M[i][j] = this->m_M[i][j] / val;
        }
    }

    return res;
}

Matrix Matrix::operator+(Matrix& rhs)
{
    if ((this->m_nRows != rhs.m_nRows) || (this->m_nCols != rhs.m_nCols))
        throw invalid_argument("Matrix + operator: dimensions must match");

    Matrix res(this->m_nRows, this->m_nCols);

    for (size_t i = 0; i < this->m_nRows; i++)
        for (size_t j = 0; j < this->m_nCols; j++)
            res.m_M[i][j] = this->m_M[i][j] + rhs.m_M[i][j];

    return res;
}

Matrix Matrix::operator+(const double value)
{
    Matrix res(this->m_nRows, this->m_nCols);

    for (size_t i = 0; i < this->m_nRows; i++)
        for (size_t j = 0; j < this->m_nCols; j++)
            res.m_M[i][j] = this->m_M[i][j] + value;

    return res;
}

Matrix Matrix::operator+(vector<double> vect)
{
    if (vect.size() != this->m_nRows)
        throw invalid_argument("Matrix + operator with vector: vector size must match the number of rows in the matrix");

    Matrix res(this->m_nRows, this->m_nCols);

    for (size_t i = 0; i < this->m_nRows; ++i) {
        for (size_t j = 0; j < this->m_nCols; ++j) {
            res.m_M[i][j] = this->m_M[i][j] + vect[i];
        }
    }
    return res;
}

Matrix Matrix::operator-(Matrix& rhs)
{
    if ((this->m_nRows != rhs.m_nRows) || (this->m_nCols != rhs.m_nCols))
        throw invalid_argument("Matrix + operator: dimensions must match");

    Matrix res(this->m_nRows, this->m_nCols);

    for (size_t i = 0; i < this->m_nRows; i++)
        for (size_t j = 0; j < this->m_nCols; j++)
            res.m_M[i][j] = this->m_M[i][j] - rhs.m_M[i][j];

    return res;
}

Matrix Matrix::operator-(const double value)
{
    Matrix res(this->m_nRows, this->m_nCols);

    for (size_t i = 0; i < this->m_nRows; i++)
        for (size_t j = 0; j < this->m_nCols; j++)
            res.m_M[i][j] = this->m_M[i][j] - value;

    return res;
}

Matrix Matrix::operator-()
{
    Matrix res(this->m_nRows, this->m_nCols);

    for (size_t i = 0; i < this->m_nRows; i++)
        for (size_t j = 0; j < this->m_nCols; j++)
            res.m_M[i][j] = -this->m_M[i][j];

    return res;
}

Matrix Matrix::operator*(Matrix& rhs)
{
    if (this->m_nCols != rhs.m_nRows)
        throw invalid_argument("Matrix * operator: incompatible dimensions");

    Matrix res(this->m_nRows, rhs.m_nCols);

    for (size_t i = 0; i < this->m_nRows; i++)
    {
        for (size_t j = 0; j < rhs.m_nCols; j++)
        {
            for (size_t k = 0; k < this->m_nCols; k++)
            {
                res.m_M[i][j] += this->m_M[i][k] * rhs.m_M[k][j];
            }
        }
    }

    return res;
}

vector<double> Matrix::operator*(vector<double>& vect)
{
    if (this->m_nRows != vect.size()) {
        throw std::invalid_argument("Matrix * operator with vector: dimensions must match");
    }

    vector<double> res(m_nRows, 0.0);

    for (size_t i = 0; i < m_nRows; ++i) {
        for (size_t j = 0; j < m_nCols; ++j) {
            res[i] += this->m_M[i][j] * vect[j];
        }
    }
    return res;
}

Matrix Matrix::operator*(const double value)
{
    Matrix res(this->m_nRows, this->m_nCols);

    for (size_t i = 0; i < this->m_nRows; i++)
    {
        for (size_t j = 0; j < this->m_nCols; j++)
        {
            res.m_M[i][j] = this->m_M[i][j] * value;
        }
    }

    return res;
}

// Accessors
double Matrix::getNCols() { return this->m_nCols; }
double Matrix::getNRows() { return this->m_nRows; }

// Mutators
void Matrix::setMatrix(vector<vector<double> >& values)
{
    if((values.size() != this->m_M.size()) || (values[0].size() != this->m_M[0].size()))
        throw invalid_argument("Dimensions must match");
    
    for(size_t i = 0; i < values.size(); i++)
    {
        for(size_t j = 0; j < values[0].size(); j++)
        {
            this->m_M[i][j] = values[i][j];
        }
    }
}

// Matrix Methods
void Matrix::display()
{
    for(size_t i = 0; i < this->m_nRows; i++)
    {
        for(size_t j = 0; j < this->m_nCols; j++)
        {
            cout << this->m_M[i][j] << " ";
        }
        cout << endl;
    }
}

double Matrix::determinant()
{
    if(this->m_nRows != this->m_nCols)
        throw invalid_argument("Matrix must be a square Matrix");
    
    size_t n = this->m_nRows;

    if(n == 1)
        return this->m_M[0][0];
    if(n == 2)
        return this->m_M[0][0] * this->m_M[1][1] - this->m_M[0][1] * this->m_M[1][0];

    double det = 0;
    Matrix subMatrix(n - 1, n - 1);

    for(size_t x = 0; x < n; x++)
    {
        size_t sub_i = 0;
        for(size_t i = 1; i < n; i++)
        {
            size_t sub_j = 0;
            for(size_t j = 0; j < n; j++)
            {
                if(j == x)
                    continue;
                
                subMatrix.m_M[sub_i][sub_j] = this->m_M[i][j];
                sub_j++;
            }
            sub_i++;
        }

        det = det + pow(-1, x) * this->m_M[0][x] * subMatrix.determinant();
    }

    return det;
}

Matrix Matrix::cofactorMatrix()
{
    size_t n = this->m_nRows;
    Matrix cofactorMat(n, n);

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            Matrix subMatrix(n - 1, n - 1);
            size_t sub_i = 0;

            for (size_t row = 0; row < n; row++)
            {
                if (row == i)
                    continue;

                size_t sub_j = 0;

                for (size_t col = 0; col < n; col++)
                {
                    if (col == j)
                        continue;

                    subMatrix.m_M[sub_i][sub_j] = this->m_M[row][col];
                    sub_j++;
                }
                sub_i++;
            }

            int sign = ((i + j) % 2 == 0) ? 1 : -1;
            cofactorMat.m_M[i][j] = sign * subMatrix.determinant();
        }
    }

    return cofactorMat;
}

Matrix Matrix::transpose()
{
    Matrix transposed(this->m_nRows, this->m_nCols);
    for(size_t i = 0; i < this->m_nRows; i++)
    {
        for(size_t j = 0; j < this->m_nCols; j++)
        {
            transposed.m_M[i][j] = this->m_M[j][i];
        }
    }
    
    return transposed;
}

Matrix Matrix::invert()
{
    double det = this->determinant();
    
    if(det == 0)
        throw invalid_argument("Matrix is not invertible : determinant is 0");

    Matrix adjugate = this->cofactorMatrix().transpose();
    return adjugate/det;
}

