
#include "../perfgaslib/perfgas.hpp"

int main() {
    
    Matrix A = {
        {0, 2, 1},
        {1, 1, 1},
        {2, 1, 0}
    };

    Vector B = {4, 6, 7};
    Vector x = B/A; 
    // displayVector(x); 

    Vector check = A * x;
    // displayVector(check); 

    return 0;
}

