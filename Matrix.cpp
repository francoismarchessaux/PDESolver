#include <cmath>
#include <iostream>
#include <vector>
#include "Matrix.h"

using namespace std;

// Default constructor
Matrix::Matrix() : m_nCols(3), m_nRows(3), m_M({{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}) { }

// Initializer Constructors
Matrix::Matrix(size_t nRows, size_t nCols) : m_nCols(nCols), m_nRows(nRows) { m_M.resize(m_nRows, vector<double>(m_nCols, 0.0)); }

Matrix::Matrix(const vector<vector<double> >& values) : m_nCols(values[0].size()), m_nRows(values.size()), m_M(values) { }

// Operators
vector<double>& Matrix::operator[](size_t i) { return m_M[i]; }

Matrix Matrix::operator/(double val) const
{
    Matrix res(m_nRows, m_nCols);

    for(size_t i = 0; i < m_nRows; i++)
        for(size_t j = 0; j < m_nCols; j++)
            res.m_M[i][j] = m_M[i][j] / val;

    return res;
}

Matrix Matrix::operator+(const Matrix& rhs) const
{
    if ((m_nRows != rhs.m_nRows) || (m_nCols != rhs.m_nCols))
        throw invalid_argument("Matrix + operator: dimensions must match");

    Matrix res(m_nRows, m_nCols);

    for (size_t i = 0; i < m_nRows; i++)
        for (size_t j = 0; j < m_nCols; j++)
            res.m_M[i][j] = m_M[i][j] + rhs.m_M[i][j];

    return res;
}

Matrix Matrix::operator+(double value) const
{
    Matrix res(m_nRows, m_nCols);

    for (size_t i = 0; i < m_nRows; i++)
        for (size_t j = 0; j < m_nCols; j++)
            res.m_M[i][j] = m_M[i][j] + value;

    return res;
}

Matrix Matrix::operator+(const vector<double> vect) const
{
    if (vect.size() != m_nRows)
        throw invalid_argument("Matrix + operator with vector: vector size must match the number of rows in the matrix");

    Matrix res(m_nRows, m_nCols);

    for (size_t i = 0; i < m_nRows; ++i)
        for (size_t j = 0; j < m_nCols; ++j)
            res.m_M[i][j] = m_M[i][j] + vect[i];

    return res;
}

Matrix Matrix::operator-(const Matrix& rhs) const
{
    if ((m_nRows != rhs.m_nRows) || (m_nCols != rhs.m_nCols))
        throw invalid_argument("Matrix + operator: dimensions must match");

    Matrix res(m_nRows, m_nCols);

    for (size_t i = 0; i < m_nRows; i++)
        for (size_t j = 0; j < m_nCols; j++)
            res.m_M[i][j] = m_M[i][j] - rhs.m_M[i][j];

    return res;
}

Matrix Matrix::operator-(double value) const
{
    Matrix res(m_nRows, m_nCols);

    for (size_t i = 0; i < m_nRows; i++)
        for (size_t j = 0; j < m_nCols; j++)
            res.m_M[i][j] = m_M[i][j] - value;

    return res;
}

Matrix Matrix::operator-() const
{
    Matrix res(m_nRows, m_nCols);

    for (size_t i = 0; i < m_nRows; i++)
        for (size_t j = 0; j < m_nCols; j++)
            res.m_M[i][j] = -m_M[i][j];

    return res;
}

Matrix Matrix::operator*(const Matrix& rhs) const
{
    if (m_nCols != rhs.m_nRows)
        throw invalid_argument("Matrix * operator: incompatible dimensions");

    Matrix res(m_nRows, rhs.m_nCols);

    for (size_t i = 0; i < m_nRows; i++)
        for (size_t j = 0; j < rhs.m_nCols; j++)
            for (size_t k = 0; k < m_nCols; k++)
                res.m_M[i][j] += m_M[i][k] * rhs.m_M[k][j];

    return res;
}

vector<double> Matrix::operator*(const vector<double>& vect) const
{
    if (m_nRows != vect.size())
        throw std::invalid_argument("Matrix * operator with vector: dimensions must match");

    vector<double> res(m_nRows, 0.0);

    for (size_t i = 0; i < m_nRows; ++i)
        for (size_t j = 0; j < m_nCols; ++j)
            res[i] += m_M[i][j] * vect[j];

    return res;
}

Matrix Matrix::operator*(double value) const
{
    Matrix res(m_nRows, m_nCols);

    for (size_t i = 0; i < m_nRows; i++)
        for (size_t j = 0; j < m_nCols; j++)
            res.m_M[i][j] = m_M[i][j] * value;

    return res;
}

// Accessors
size_t Matrix::getNCols() { return m_nCols; }
size_t Matrix::getNRows() { return m_nRows; }

// Mutators
void Matrix::setMatrix(const vector<vector<double> >& values)
{
    if((values.size() != m_M.size()) || (values[0].size() != m_M[0].size()))
        throw invalid_argument("Dimensions must match");
    
    for(size_t i = 0; i < values.size(); i++)
        for(size_t j = 0; j < values[0].size(); j++)
            m_M[i][j] = values[i][j];
}

// Matrix Methods
void Matrix::display() const
{
    for(size_t i = 0; i < m_nRows; i++)
    {
        for(size_t j = 0; j < m_nCols; j++)
            cout << m_M[i][j] << " ";
        cout << endl;
    }
}

Matrix Matrix::invert() const
{
    if (m_nRows != m_nCols)
        throw std::invalid_argument("Matrix must be square for inversion");

    size_t n = m_nRows;
    Matrix result(n, n);
    Matrix augmentedMatrix(n, 2 * n);

    // Initialize augmented matrix
    for (size_t i = 0; i < n; ++i) 
    {
        for (size_t j = 0; j < n; ++j) 
        {
            augmentedMatrix[i][j] = m_M[i][j];
            if(i == j)
                augmentedMatrix[i][j + n] = 1.0;
            else
                augmentedMatrix[i][j + n] = 0.0;
        }
    }

    // Gaussian elimination
    for (size_t i = 0; i < n; ++i) 
    {
        double pivot = augmentedMatrix[i][i];

        if (pivot == 0.0)
            throw std::invalid_argument("Matrix is not invertible");

        // Scale the current row to make the pivot 1
        for (size_t j = 0; j < 2 * n; ++j)
            augmentedMatrix[i][j] /= pivot;

        // Eliminate other rows
        for (size_t k = 0; k < n; ++k) 
        {
            if (k != i) 
            {
                double factor = augmentedMatrix[k][i];

                for (size_t j = 0; j < 2 * n; ++j)
                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
            }
        }
    }

    // Extract the inverted matrix from the augmented matrix
    for (size_t i = 0; i < n; ++i) 
        for (size_t j = 0; j < n; ++j) 
            result[i][j] = augmentedMatrix[i][j + n];

    return result;
}
