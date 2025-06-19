
#include "../perfgaslib/perfgas.hpp"

int main() {
    
    const int Nx = 200, Ny = 100;

    RampGrid grid(Nx, Ny, 3, 1, 20);

    string filename = "../plotfiles/gridplot.dat";

    output_grid(grid, Nx, Ny, filename);

    return 0;
}
