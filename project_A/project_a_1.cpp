#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;

void write_xyz(vector <vector <double>> coords, string atom="Si", string filename="cell.xyz", string comment="\n") {
    ofstream file;
    file.open(filename);
    file << coords.size() << "\n";
    file << comment;
    for (auto xyz: coords) {
        file << atom << "\t" << xyz[0] << "\t" << xyz[1] << "\t" << xyz[2] << "\n";
    }
    file.close();
    return;
}

vector <vector <double>> get_simplecubic(int nx, int ny, int nz, double lattice_cnst) {
    vector <vector <double>> coords;
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                coords.push_back({i*lattice_cnst, j*lattice_cnst, k*lattice_cnst});
            }
        }
    }
    return coords;
}

vector <vector <double>> get_bcc(int nx, int ny, int nz, double lattice_cnst) {
    vector <vector <double>> coords;
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                coords.push_back({i*lattice_cnst, j*lattice_cnst, k*lattice_cnst});
                coords.push_back({(i + 0.5)*lattice_cnst, (j + 0.5)*lattice_cnst, (k + 0.5)*lattice_cnst});
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
                coords.push_back({i*lattice_cnst, j*lattice_cnst, k*lattice_cnst});
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
                coords.push_back({i*lattice_cnst, j*lattice_cnst, k*lattice_cnst});
                coords.push_back({(i + 0.5)*lattice_cnst, (j + 0.5)*lattice_cnst, k*lattice_cnst});
                coords.push_back({(i + 0.5)*lattice_cnst, j*lattice_cnst, (k + 0.5)*lattice_cnst});
                coords.push_back({i*lattice_cnst, (j + 0.5)*lattice_cnst, (k + 0.5)*lattice_cnst});
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
    int nx, ny, nz;
    double a;
    cout << "Please enter Nx Ny Nz (for cell dimension Nx by Ny by Nz): ";
    cin >> nx >> ny >> nz;
    cout << "Please enter lattice constant (float): ";
    cin >> a;
    int num_atoms = 0;
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