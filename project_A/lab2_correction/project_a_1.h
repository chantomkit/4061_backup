#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;

void write_xyz(vector <vector <double> > coords, string atom="Si", string filename="cell.xyz", string comment="\n") {
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

vector <vector <double> > get_simplecubic(int nx, int ny, int nz, double lattice_cnst, double shiftx = 0, double shifty = 0, double shiftz = 0) {
    // coords is a 2D array for storing (x,y,z) of each atom at each row
    vector <vector <double> > coords;
    // iterates through all x, y, z directions with given periods
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                // looping through the lattice points and scaled by the lattice constant
                coords.push_back({i*lattice_cnst+shiftx, j*lattice_cnst+shifty, k*lattice_cnst+shiftz}); // this adds cubic point (corner)
            }
        }
    }
    return coords;
}

vector <vector <double> > get_bcc(int nx, int ny, int nz, double lattice_cnst, double shiftx = 0, double shifty = 0, double shiftz = 0) {
    // coords is a 2D array for storing (x,y,z) of each atom at each row
    vector <vector <double> > coords;
    // iterates through all x, y, z directions with given periods
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                // looping through the lattice points and scaled by the lattice constant
                coords.push_back({i*lattice_cnst+shiftx, j*lattice_cnst+shifty, k*lattice_cnst+shiftz}); // this adds cubic point (corner)
                coords.push_back({(i + 0.5)*lattice_cnst+shiftx, (j + 0.5)*lattice_cnst+shifty, (k + 0.5)*lattice_cnst+shiftz}); // this adds the body centered point
            }
        }
    }
    return coords;
}

vector <vector <double> > get_fcc(int nx, int ny, int nz, double lattice_cnst, double shiftx = 0, double shifty = 0, double shiftz = 0) {
    vector <vector <double> > coords;
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                // looping through the lattice points and scaled by the lattice constant
                coords.push_back({i*lattice_cnst+shiftx, j*lattice_cnst+shifty, k*lattice_cnst+shiftz}); // this adds cubic point (corner)
                // These add the face centered points (3 points only), such that total atoms in a cell = 4
                coords.push_back({(i + 0.5)*lattice_cnst+shiftx, (j + 0.5)*lattice_cnst+shifty, k*lattice_cnst+shiftz});
                coords.push_back({(i + 0.5)*lattice_cnst+shiftx, j*lattice_cnst+shifty, (k + 0.5)*lattice_cnst+shiftz});
                coords.push_back({i*lattice_cnst+shiftx, (j + 0.5)*lattice_cnst+shifty, (k + 0.5)*lattice_cnst+shiftz});
            }
        }
    }
    return coords;
}

vector <vector <double> > get_diamond(int nx, int ny, int nz, double lattice_cnst, double shiftx = 0, double shifty = 0, double shiftz = 0) {
    vector <vector <double> > coords;
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                // looping through the lattice points and scaled by the lattice constant
                coords.push_back({i*lattice_cnst+shiftx, j*lattice_cnst+shifty, k*lattice_cnst}); // this adds cubic point (corner)
                // These add face center points (3 points only)
                coords.push_back({(i + 0.5)*lattice_cnst+shiftx, (j + 0.5)*lattice_cnst+shifty, k*lattice_cnst+shiftz});
                coords.push_back({(i + 0.5)*lattice_cnst+shiftx, j*lattice_cnst+shifty, (k + 0.5)*lattice_cnst+shiftz});
                coords.push_back({i*lattice_cnst+shiftx, (j + 0.5)*lattice_cnst+shifty, (k + 0.5)*lattice_cnst+shiftz});
                // These add diamond structure points (4 points only), such that total 8 atoms in a cell
                coords.push_back({(i + 0.25)*lattice_cnst+shiftx, (j + 0.25)*lattice_cnst+shifty, (k + 0.25)*lattice_cnst+shiftz});
                coords.push_back({(i + 0.75)*lattice_cnst+shiftx, (j + 0.75)*lattice_cnst+shifty, (k + 0.25)*lattice_cnst+shiftz});
                coords.push_back({(i + 0.25)*lattice_cnst+shiftx, (j + 0.75)*lattice_cnst+shifty, (k + 0.75)*lattice_cnst+shiftz});
                coords.push_back({(i + 0.75)*lattice_cnst+shiftx, (j + 0.25)*lattice_cnst+shifty, (k + 0.75)*lattice_cnst+shiftz});
            }
        }
    }
    return coords;
}