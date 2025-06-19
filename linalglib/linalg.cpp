// Linear Algebra Library.cpp : Defines the functions for the static library.
//

#include "linalg.hpp" 


// Creates an identity matrix of given size
Matrix identity(int a) {
    Matrix I(a, Vector(a, 0.0));
    for (int i = 0; i < 4; ++i) {
        I[i][i] = 1.0;
    }
    return I;
}

// This function creates a vector full of random variables.
Vector random(int a) {
    Vector result(a);

    for (int j = 0; j < 4; ++j) {
        result[j] = rand();
    }
    return result;
}

// This function creates a matrix full of random variables.
Matrix random(int a, int b) {

    Matrix result(a, Vector(b, 0.0));

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result[i][j] = rand();
        }
    }

    return result;
}

void displayVector(const Vector& A) {
    for (int j = 0; j < A.size(); ++j) {
        cout << A[j] << " ";
    }
}

void displayMatrix(const Matrix& A) {
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[0].size(); ++j) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}

Matrix outerProduct(const Vector& a, const Vector& b) {
    Matrix result(a.size(), Vector(a.size()));

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result[i][j] = a[i] * b[j];
        }
    }
    return result;
}

Vector operator+(const Vector& v1, const Vector& v2) {

    Vector result = zeros(v1.size());
    for (int i = 0; i < 4; i++) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

double operator*(const Vector& A, const Vector& B) {
    double result = 0.0;
    for (int i = 0; i < 4; i++) {
        result += A[i] * B[i];
    }
    return result;
}

Vector operator-(const Vector& v1, const Vector& v2) {

    Vector result(v1.size());
    for (size_t i = 0; i < 4; i++) {
        result[i] = v1[i] - v2[i];
    }
    return result;
}

Vector operator*(const Vector& A, const double& s) {

    Vector result(A.size());

    for (int i = 0; i < 4; ++i) {
        result[i] = s * A[i]; // Perform scalar multiplication
    }

    return result;
}

Vector operator*(const double& s, const Vector& A) {
    return A * s;
}

Vector operator/(const Vector& v1, const double& s) {

    Vector result(v1.size());
    for (int i = 0; i < 4; i++) {
        result[i] = v1[i] / s;
    }

    return result;
}

Vector operator*(const Matrix& A, const Vector& B) {

    Vector result = zeros(B.size());

    for (int i = 0; i < 4; ++i) {
        for (int k = 0; k < 4; ++k) {
            result[i] += A[i][k] * B[k];
        }
    }
    return result;
}

Matrix operator/(const Matrix& A, const double& s) {

    Matrix result(A.size(), Vector(A[0].size()));

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result[i][j] = A[i][j] / s;
        }
    }

    return result;
}

Matrix operator+(const Matrix& A, const Matrix& B) {

    Matrix result(A.size(), Vector(A[0].size()));

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }

    return result;
}

Matrix operator-(const Matrix& A, const Matrix& B) {

    // Initialize result matrix with proper size
    Matrix result(A.size(), Vector(A[0].size()));

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }

    return result;
}

Matrix operator*(const double& s, const Matrix& A) {

    Matrix result(A.size(), Vector(A[0].size()));

    // Multiply each element of the matrix by the scalar
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result[i][j] = s * A[i][j]; // Perform scalar multiplication
        }
    }

    return result;
}

Matrix operator*(const Matrix& A, const double& s) {
    return s * A;
}

Matrix operator*(const Matrix& A, const Matrix& B) {

    Matrix result(A.size(), Vector(B[0].size()));

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}


// General LU decomposition with permutation matrix for partial pivoting //
void LUDecomposition(const Matrix& A, Matrix& L, Matrix& U, vector<int>& P) {
    int n = A.size();
    L = zeros(n, n);
    U = A;
    P.resize(n);
    for (int i = 0; i < n; i++) P[i] = i;


    for (int i = 0; i < n; i++) {
        // Pivot
        double maxVal = std::abs(U[i][i]);
        int pivot = i;
        for (int j = i + 1; j < n; j++) {
            if (std::abs(U[j][i]) > maxVal) {
                maxVal = std::abs(U[j][i]);
                pivot = j;
            }
        }


        if (pivot != i) {
            std::swap(U[i], U[pivot]);
            std::swap(P[i], P[pivot]);
            for (int k = 0; k < i; k++) {
                std::swap(L[i][k], L[pivot][k]);
            }
        }


        L[i][i] = 1.0;
        for (int j = i + 1; j < n; j++) {
            if (U[i][i] == 0.0)
                throw std::runtime_error("Singular matrix detected!");


            double factor = U[j][i] / U[i][i];
            L[j][i] = factor;
            for (int k = i; k < n; k++) {
                U[j][k] -= factor * U[i][k];
            }
        }
    }
}

// Matrix - Matrix "division" //
Matrix applyPermutation(const Matrix& B, const vector<int>& P) {
    int n = P.size();
    int m = B[0].size();
    Matrix Pb(n, Vector(m));
    for (int i = 0; i < n; i++) {
        Pb[i] = B[P[i]];
    }
    return Pb;
}
Matrix forwardSubstitution(const Matrix& L, const Matrix& B) {
    int n = L.size();
    int m = B[0].size(); // Number of right-hand sides
    Matrix Y = zeros(n, m);


    for (int i = 0; i < n; i++) {
        for (int col = 0; col < m; col++) {
            Y[i][col] = B[i][col];
            for (int j = 0; j < i; j++) {
                Y[i][col] -= L[i][j] * Y[j][col];
            }
        }
    }
    return Y;
}
Matrix backwardSubstitution(const Matrix& U, const Matrix& Y) {
    int n = U.size();
    int m = Y[0].size();
    Matrix X = zeros(n, m);


    for (int i = n - 1; i >= 0; i--) {
        for (int col = 0; col < m; col++) {
            if (U[i][i] == 0.0) {
                throw std::runtime_error("Singular matrix detected in back-substitution.");
            }
            X[i][col] = Y[i][col];
            for (int j = i + 1; j < n; j++) {
                X[i][col] -= U[i][j] * X[j][col];
            }
            X[i][col] /= U[i][i];
        }
    }
    return X;
}
Matrix operator/(const Matrix& B, const Matrix& A) {


    Matrix L, U;
    std::vector<int> P;


    LUDecomposition(A, L, U, P);          // Pivoting version
    Matrix Pb = applyPermutation(B, P);   // Permute right-hand sides
    Matrix Y = forwardSubstitution(L, Pb);
    return backwardSubstitution(U, Y);
}


// Matrix - Vector "division" //
Vector forwardSubstitution(const Matrix& L, const Vector& B) {
    int n = L.size();
    Vector Y(n);
    for (int i = 0; i < n; i++) {
        Y[i] = B[i];
        for (int j = 0; j < i; j++) {
            Y[i] -= L[i][j] * Y[j];
        }
    }
    return Y;
}
Vector backwardSubstitution(const Matrix& U, const Vector& Y) {
    int n = U.size();
    Vector X(n);
    for (int i = n - 1; i >= 0; i--) {
        if (U[i][i] == 0.0)
            throw std::runtime_error("Singular matrix in back-substitution.");
        X[i] = Y[i];
        for (int j = i + 1; j < n; j++) {
            X[i] -= U[i][j] * X[j];
        }
        X[i] /= U[i][i];
    }
    return X;
}
Vector applyPermutation(const Vector& B, const vector<int>& P) {
    int n = B.size();
    Vector out(n);
    for (int i = 0; i < n; i++)
        out[i] = B[P[i]];
    return out;
}
Vector operator/(const Vector& B, const Matrix& A) {
    Matrix L, U;
    std::vector<int> P;
    LUDecomposition(A, L, U, P);
    Vector Pb = applyPermutation(B, P);
    Vector Y = forwardSubstitution(L, Pb);
    return backwardSubstitution(U, Y);
}
