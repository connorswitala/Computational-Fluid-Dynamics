
#include "../2D_Perf_Gas_FVS/perf_gas_fvs.hpp"

int main() {
    
    const int Nx = 200, Ny = 100;  
    RampGrid grid(Nx, Ny, 1, 1, 1, 0.65, 15);  

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            if (isnan(grid.jNorms(i,j).x)) {
                cout << "ISNAN" << endl; 
            }
            if (isinf(grid.jNorms(i,j).y)) { 
                cout << "ISINF" << endl; 
            }
        }
    }

    return 0;
}

