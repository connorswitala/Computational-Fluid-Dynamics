#pragma once

#include <iostream>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <stdexcept>  // For exception handling
#include <chrono>
#include <unordered_map>
#include <string>
#include <array>

using namespace std;
constexpr double pi = 3.141592653;



typedef vector<double> Vector;
typedef vector<Vector> Matrix;
typedef vector<Matrix> Tensor;
typedef vector<Tensor> Tesseract;


inline Vector zeros(int a) {
    return Vector(a, 0.0);
}

inline Matrix zeros(int a, int b) {
    return Matrix(a, Vector(b, 0.0));
}

inline Vector ones(int a) {
    return Vector{ 1, 1, 1, 1 };
}

inline Matrix ones(int a, int b) {
    return Matrix(a, Vector(b, 1.0));
}

Matrix identity(int a);
Vector random(int a);
Matrix random(int a, int b);

void displayVector(const Vector& A);
void displayMatrix(const Matrix& A);

double operator*(const Vector& A, const Vector& B);
Vector operator+(const Vector& v1, const Vector& v2);
Vector operator-(const Vector& v1, const Vector& v2);
Vector operator*(const Vector& A, const double& s);
Vector operator*(const double& s, const Vector& A);
Vector operator/(const Vector& v1, const double& s);


Matrix outerProduct(const Vector& a, const Vector& b);

Vector operator*(const Matrix& A, const Vector& B);


Matrix operator+(const Matrix& A, const Matrix& B);
Matrix operator-(const Matrix& A, const Matrix& B);
Matrix operator*(const Matrix& A, const double& s);
Matrix operator*(const double& s, const Matrix& A);
Matrix operator*(const Matrix& A, const Matrix& B);
Matrix operator/(const Matrix& A, const double& s);


void LUDecomposition(const Matrix& A, Matrix& L, Matrix& U);
Vector forwardSubstitution(const Matrix& L, const Vector& B);
Vector backwardSubstitution(const Matrix& U, const Vector& Y);
Vector operator/(const Vector& B, const Matrix& A);


Matrix forwardSubstitution(const Matrix& L, const Matrix& b);
Matrix backwardSubstitution(const Matrix& U, const Matrix& y);
Matrix operator/(const Matrix& b, const Matrix& A);