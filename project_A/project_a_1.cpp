#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;

void write_xyz(vector <vector <double>> coords, string atom="Si", string filename="cell.xyz", string comment="\n") {
    ofstream file;
    file.open(filename);
    file << coords.size() << "\n"; // get length of coords array as the total number of atoms of the crystal
    file << comment; // comment line, left empty by default
    for (auto xyz: coords) {
        file << atom << "\t" << xyz[0] << "\t" << xyz[1] << "\t" << xyz[2] << "\n"; // line by line writing coords of each atom
    }
    file.close();
    return;
}

vector <vector <double>> get_simplecubic(int nx, int ny, int nz, double lattice_cnst) {
    // coords is a 2D array for storing (x,y,z) of each atom at each row
    vector <vector <double>> coords;
    // iterates through all x, y, z directions with given periods
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                // looping through the lattice points and scaled by the lattice constant
                coords.push_back({i*lattice_cnst, j*lattice_cnst, k*lattice_cnst}); // this adds cubic point (corner)
            }
        }
    }
    return coords;
}

vector <vector <double>> get_bcc(int nx, int ny, int nz, double lattice_cnst) {
    // coords is a 2D array for storing (x,y,z) of each atom at each row
    vector <vector <double>> coords;
    // iterates through all x, y, z directions with given periods
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                // looping through the lattice points and scaled by the lattice constant
                coords.push_back({i*lattice_cnst, j*lattice_cnst, k*lattice_cnst}); // this adds cubic point (corner)
                coords.push_back({(i + 0.5)*lattice_cnst, (j + 0.5)*lattice_cnst, (k + 0.5)*lattice_cnst}); // this adds the body centered point
            }
        }
    }
    return coords;
}

vector <vector <double>> get_fcc(int nx, int ny, int nz, double lattice_cnst) {
    vector <vector <double>> coords;
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                // looping through the lattice points and scaled by the lattice constant
                coords.push_back({i*lattice_cnst, j*lattice_cnst, k*lattice_cnst}); // this adds cubic point (corner)
                // These add the face centered points (3 points only), such that total atoms in a cell = 4
                coords.push_back({(i + 0.5)*lattice_cnst, (j + 0.5)*lattice_cnst, k*lattice_cnst});
                coords.push_back({(i + 0.5)*lattice_cnst, j*lattice_cnst, (k + 0.5)*lattice_cnst});
                coords.push_back({i*lattice_cnst, (j + 0.5)*lattice_cnst, (k + 0.5)*lattice_cnst});
            }
        }
    }
    return coords;
}

vector <vector <double>> get_diamond(int nx, int ny, int nz, double lattice_cnst) {
    vector <vector <double>> coords;
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                // looping through the lattice points and scaled by the lattice constant
                coords.push_back({i*lattice_cnst, j*lattice_cnst, k*lattice_cnst}); // this adds cubic point (corner)
                // These add face center points (3 points only)
                coords.push_back({(i + 0.5)*lattice_cnst, (j + 0.5)*lattice_cnst, k*lattice_cnst});
                coords.push_back({(i + 0.5)*lattice_cnst, j*lattice_cnst, (k + 0.5)*lattice_cnst});
                coords.push_back({i*lattice_cnst, (j + 0.5)*lattice_cnst, (k + 0.5)*lattice_cnst});
                // These add diamond structure points (4 points only), such that total 8 atoms in a cell
                coords.push_back({(i + 0.25)*lattice_cnst, (j + 0.25)*lattice_cnst, (k + 0.25)*lattice_cnst});
                coords.push_back({(i + 0.75)*lattice_cnst, (j + 0.75)*lattice_cnst, (k + 0.25)*lattice_cnst});
                coords.push_back({(i + 0.25)*lattice_cnst, (j + 0.75)*lattice_cnst, (k + 0.75)*lattice_cnst});
                coords.push_back({(i + 0.75)*lattice_cnst, (j + 0.25)*lattice_cnst, (k + 0.75)*lattice_cnst});
            }
        }
    }
    return coords;
}

int main() {
    int nx, ny, nz; // x,y,z dimensions
    double a; // lattice constant
    cout << "Please enter Nx Ny Nz (for cell dimension Nx by Ny by Nz): ";
    cin >> nx >> ny >> nz;
    cout << "Please enter lattice constant (float): ";
    cin >> a;
    cout << nx << " " << ny << " " << nz << " " << a << endl;
    vector <vector <double>> simple_cubic = get_simplecubic(nx, ny, nz, a);
    write_xyz(simple_cubic, "Si", "sc.xyz");
    vector <vector <double>> bcc = get_bcc(nx, ny, nz, a);
    write_xyz(bcc, "Si", "bcc.xyz");
    vector <vector <double>> fcc = get_fcc(nx, ny, nz, a);
    write_xyz(fcc, "Si", "fcc.xyz");
    vector <vector <double>> diamond = get_diamond(nx, ny, nz, a);
    write_xyz(diamond, "Si", "diamond.xyz");
    return 0;
}