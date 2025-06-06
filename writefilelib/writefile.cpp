#include "writefile.hpp"

Vector convert(const Vector& U) {
    double gam = 1.4; 

	static Vector V(4);
	V[0] = U[0];
	V[1] = U[1] / U[0];
	V[2] = U[2] / U[0];
	V[3] = (U[3] - 0.5 * V[0] * (V[1] * V[1] + V[2] * V[2])) * (gam - 1);

	return V;
} 
Vector convert_real(const Vector& U, ThermoEntry thermo) {

	static Vector V(4);
	V[0] = U[0];
	V[1] = U[1] / U[0];
	V[2] = U[2] / U[0];
	V[3] = (U[3] - 0.5 * V[0] * (V[1] * V[1] + V[2] * V[2])) * (thermo.gamma - 1); 

	return V;
} 

void output_grid(Grid& grid, const int& Nx, const int& Ny, string& filename) {
    

    ofstream file(filename);
    file << "Variables = \" x points\", \"y points\" \n";    
    file << "ZONE T=\"Grid\", I=" << Nx << ", J=" << Ny << ", F=POINT\n";



    for (int i = 0; i <= Nx; ++i) {
        for (int j = 0; j <= Ny; ++j) {

            file << grid.Vertex(i,j).x << " " << grid.Vertex(i,j).y << endl;             
        }
    }

    file.close(); 
}

void cfd_contour(const int Nx, const int Ny, Tensor& U, Grid& grid, string filename) {

    double gam = 1.4; 
    double R = 287.0; 

    ofstream file(filename);
    file << "Variables = \" x points\", \"y points\", \"density\", \"pressure\", \"temperature\", \"mach\", \"u velocity\", \"v velocity\" \n";    
    file << "ZONE T=\"Flow Field\", I=" << Nx << ", J=" << Ny << ", F=POINT\n";


    Vector primitives(4); 
    double temp, mach, a;

    for (int i = 1; i < Nx + 1; ++i) {
        for (int j = 1; j < Ny + 1; ++j) {
            primitives = convert(U[i][j]);  
            temp = primitives[3]/(R * primitives[0]);
            a = sqrt(gam * R * temp);
            mach = sqrt(primitives[1] * primitives[1] + primitives[2] * primitives[2]) / a;  

            file << grid.Center(i-1,j-1).x << " " << grid.Center(i-1,j-1).y << " " << primitives[0] 
            << " " << primitives[3] << " " << temp << " " << mach << " " << primitives[1] << " " <<
            primitives[2] << endl; 
            
        }
    }

    file.close(); 
}

void cfd_real_contour(const int Nx, const int Ny, Tensor& U, Grid& grid, vector<vector<ThermoEntry>>& cell_thermo, string filename) {

    ofstream file(filename);
    file << "Variables = \" x points\", \"y points\", \"density\", \"pressure\", \"temperature\", \"mach\", \"u velocity\", \"v velocity\" \n";    
    file << "ZONE T=\"Flow Field\", I=" << Nx << ", J=" << Ny << ", F=POINT\n";


    Vector primitives(4); 
    double temp, mach, a;

    for (int i = 1; i < Nx + 1; ++i) {
        for (int j = 1; j < Ny + 1; ++j) {
            primitives = convert_real(U[i][j], cell_thermo[i][j]);   
            temp = cell_thermo[i][j].p/(cell_thermo[i][j].R * cell_thermo[i][j].rho);  
            a = sqrt(cell_thermo[i][j].gamma * cell_thermo[i][j].R * cell_thermo[i][j].T);   
            mach = sqrt(primitives[1] * primitives[1] + primitives[2] * primitives[2]) / a;  

            file << grid.Center(i-1,j-1).x << " " << grid.Center(i-1,j-1).y << " " << primitives[0] 
            << " " << primitives[3] << " " << temp << " " << mach << " " << primitives[1] << " " <<
            primitives[2] << endl; 
            
        }
    }

    file.close(); 
}


void cfd_centerline(const int Nx, const int Ny, Tensor& U, Grid& grid, string filename) {

    double gam = 1.4; 
    double R = 287.0; 

    ofstream file(filename);
    file << "Variables = \"x points\", \"y points\", \"density\", \"pressure\", \"temperature\", \"mach\", \"u velocity\", \"v velocity\" \n";
    file << "Zone I=" << Nx * Ny << ", F=POINT\n";

    Vector primitives(4); 
    double temp, mach, a;
    int center = round(static_cast<double>(Ny) / 2); 

    for (int i = 0; i < Nx; ++i) {

            primitives = convert(U[i][center]); 
            temp = primitives[3]/(R * primitives[0]);
            a = sqrt(gam * R * temp);
            mach = sqrt(primitives[1]*primitives[1] + primitives[2]*primitives[2])/a; 

            file << grid.Center(i,center).x << " " << grid.Center(i,center).y << " " << primitives[0] 
            << " " << primitives[3] << " " << temp << " " << mach << " " << primitives[1] << " " <<
            primitives[2] << endl; 
    }
    

    file.close(); 
}