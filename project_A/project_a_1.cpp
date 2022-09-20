#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include "project_a_1.h"
using namespace std;

int main() {
    int nx, ny, nz; // x,y,z dimensions
    double a; // lattice constant
    cout << "Please enter Nx Ny Nz (for cell dimension Nx by Ny by Nz): ";
    cin >> nx >> ny >> nz;
    cout << "Please enter lattice constant (float): ";
    cin >> a;
    cout << nx << " " << ny << " " << nz << " " << a << endl;
    vector <vector <double> > simple_cubic = get_simplecubic(nx, ny, nz, a);
    write_xyz(simple_cubic, "Si", "sc.xyz");
    vector <vector <double> > bcc = get_bcc(nx, ny, nz, a);
    write_xyz(bcc, "Si", "bcc.xyz");
    vector <vector <double> > fcc = get_fcc(nx, ny, nz, a);
    write_xyz(fcc, "Si", "fcc.xyz");
    vector <vector <double> > diamond = get_diamond(nx, ny, nz, a);
    write_xyz(diamond, "Si", "diamond.xyz");
    return 0;
}