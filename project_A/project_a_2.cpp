#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
// Project A part 1 functions are reused
// simple cubic, bcc, fcc, diamond structure generation
#include "project_a_1.h"
using namespace std;

// Custom struct for storing neighbour records
struct neighbour_record {
    int label1;
    int label2;
    double distance;
};

// Cross product function (3-D)
vector <double> cross_prod(vector <double> u, vector <double> v) {
    double x = u[1]*v[2] - u[2]*v[1];
    double y = -u[0]*v[2] + u[2]*v[0];
    double z = u[0]*v[1] - u[1]*v[0];
    return {x,y,z};
}

// Dot product function (N-D)
double dot_prod(vector <double> u, vector <double> v) {
    double sum = 0;
    for (int i = 0; i < u.size(); i++)
    {
        sum += u[i]*v[i];
    }
    return sum;
}

// Volume (3-D)
double volume(vector <vector <double> > lattice) {
    return dot_prod(lattice[0], cross_prod(lattice[1], lattice[2]));
}

// Reciprocal vector calculation
// input 3 lattice vectors, output 3 reciprocal lattice vectors
vector <vector <double> > reciprocal_vec(vector <vector <double> > a) {
    double volume = dot_prod(a[0], cross_prod(a[1], a[2]));
    
    vector <double> b1;
    for (auto coef: cross_prod(a[1], a[2]))
    {
        b1.push_back(coef / volume);
    }
    vector <double> b2;
    for (auto coef: cross_prod(a[2], a[0]))
    {
        b2.push_back(coef / volume);
    }
    vector <double> b3;
    for (auto coef: cross_prod(a[0], a[1]))
    {
        b3.push_back(coef / volume);
    }
    return {b1, b2, b3};
}

// misc function for output
void print_1dvector(vector <double> vec) {
    for (auto x: vec)
        {
            cout << x << " ";
        }
        cout << endl;
    return;
}

// misc function for output
void print_2dvector(vector <vector <double> > vec) {
    for (auto v: vec)
    {
        print_1dvector(v);
    }
    return;
}

// PBC on fractional coords
// input a 1-D vector, output a 1-D vector with coefficients in [-0.5, 0.5]
vector <double> pbc(vector <double> d) {
    double temp;
    vector <double> res;
    for (auto x: d)
    {
        temp = fmod(x, 1.0);
        if (temp > 0.5)
        {
            temp -= 1;
        }
        else if (temp < -0.5)
        {
            temp += 1;
        }
        res.push_back(temp);
    }
    return res;
}

// True distance
// input: distance vector (1-D), real lattice vectors, reciprocal lattice vectors, print option
// output: PBC reciprocal distance vector and PBC real space distance vector
pair<vector <double>, vector <double>>true_distance_vec(vector <double> d, vector <vector <double> > real_lattice, vector <vector <double> > recip_lattice, bool print=false) {
    vector <double> recip_d, true_d;
    if (print)
        cout << "Reciprocal distance vector" << endl;
    for (auto recip_vec: recip_lattice)
    {
        // Step 1: Compute reciprocal distance vector
        // Each reciprocal lattice vector dot distance vector (1D)
        // Yields reciprocal projection of distance vector
        recip_d.push_back(dot_prod(recip_vec, d));
        if (print)
            cout << recip_d.back() << " ";
    }
    // Step 2: Apply PBC on reciprocal distance vector
    // i.e. PBC in reciprocal space
    recip_d = pbc(recip_d);
    if (print) {
        cout << endl;
        cout << "PBC reciprocal distance vector" << endl;
        print_1dvector(recip_d);
        cout << "PBC real distance vector" << endl;
    }
    for (auto real_vec: real_lattice)
    {
        // Step 3: Convert the PBC reciprocal distance vector to real space
        // Each real lattice vector dot PBC reciprocal distance vector
        // Yields PBC real space distance vector
        true_d.push_back(dot_prod(real_vec, recip_d));
        if (print)
            cout << true_d.back() << " ";
    }
    if (print)
        cout << endl;
    return make_pair(recip_d, true_d);
}

// Neighbout list function
// input: real lattice vectors, xyz coords of all atoms in the structure, distance cutoff
// output: neighbour list (Array of neighbour record struct)
vector <neighbour_record> neighbour_list(vector <vector <double> > lattice, vector <vector <double> > xyz, double dist_cutoff) {
    vector <neighbour_record> nlist;
    vector <double> dist_vec_temp, true_dist_vec_temp;
    double true_dist_temp;
    // Iterate all possible pair of atoms (without self-self counting and double counting)
    for (int i = 0; i < xyz.size(); i++)
    {
        for (int j = i; j < xyz.size(); j++)
        {
            if (i == j)
                continue;
            // Step 1: Compute distance vector of chosen pair of atoms
            dist_vec_temp = {xyz[j][0] - xyz[i][0], xyz[j][1] - xyz[i][1], xyz[j][2] - xyz[i][2]};
            // Step 2: Get true distance vector in real space (PBC applied)
            true_dist_vec_temp = true_distance_vec(dist_vec_temp, lattice, reciprocal_vec(lattice)).second;
            // Step 3: Self dot product of true distance vector to get distance squared
            true_dist_temp = dot_prod(true_dist_vec_temp, true_dist_vec_temp);
            // Step 4: if distance squared is smaller than distance cutoff squared, add the pair to neighbour list
            if (true_dist_temp < (dist_cutoff * dist_cutoff))
            {
                nlist.push_back({i, j, true_dist_temp});
            }
        }
    }
    return nlist;
}

// Write the neighbour list into a csv file
void write_neighbour_list(vector <neighbour_record> nlist, string filename="neighbour_list.csv") {
    ofstream file;
    file.open(filename);
    file << "label1,label2,distance\n";
    for (auto record: nlist) {
        file << record.label1 << "," << record.label2 << "," << record.distance << "\n";
    }
    file.close();
    return;
}

int main() {
    double lat_cnst = 1;
    vector <vector <double> > sc_vec = {{lat_cnst,0,0}, {0,lat_cnst,0}, {0,0,lat_cnst}};
    vector <vector <double> > sc_recip = reciprocal_vec(sc_vec);
    cout << "RECIPROCAL FUNCTION TESTING PART" << endl;
    cout << "Simple Cubic Primitive Vector" << endl;
    print_2dvector(sc_vec);
    cout << "Simple Cubic Reciprocal Vector" << endl;
    print_2dvector(sc_recip);
    cout << "Primitive Simple Cube Cell Volume: " << volume(sc_vec) << endl;
    cout << "Reciprocal Simple Cube Cell Volume: " << volume(sc_recip) << endl << endl;

    vector <vector <double> > bcc_vec = {{lat_cnst/2,lat_cnst/2,-lat_cnst/2}, {-lat_cnst/2,lat_cnst/2,lat_cnst/2}, {lat_cnst/2,-lat_cnst/2,lat_cnst/2}};
    vector <vector <double> > bcc_recip = reciprocal_vec(bcc_vec);
    cout << "BCC Primitive Vector" << endl;
    print_2dvector(bcc_vec);
    cout << "BCC Reciprocal Vector" << endl;
    print_2dvector(bcc_recip);
    cout << "Primitive BCC Volume: " << volume(bcc_vec) << endl;
    cout << "Reciprocal BCC Volume: " << volume(bcc_recip) << endl << endl;

    vector <vector <double> > fcc_vec = {{lat_cnst/2,lat_cnst/2,0}, {0,lat_cnst/2,lat_cnst/2}, {lat_cnst/2,0,lat_cnst/2}};
    vector <vector <double> > fcc_recip = reciprocal_vec(fcc_vec);
    cout << "FCC Primitive Vector" << endl;
    print_2dvector(fcc_vec);
    cout << "FCC Reciprocal Vector" << endl;
    print_2dvector(fcc_recip);
    cout << "Primitive FCC Volume: " << volume(fcc_vec) << endl;
    cout << "Reciprocal FCC Volume: " << volume(fcc_recip) << endl << endl;

    vector <double> dist = {1.62, -8.24, 5.36};
    vector <vector <double> > test_ucell_vec = {{2, 0, 0}, {0, 3, 0}, {0, 0, 4}};
    cout << "TRUE DISTANCE FUNCTION TESTING PART" << endl;
    cout << "Unit Cell Vector" << endl;
    print_2dvector(test_ucell_vec);
    cout << "Original Distance Vector" << endl;
    print_1dvector(dist);
    true_distance_vec(dist, test_ucell_vec, reciprocal_vec(test_ucell_vec), true);
    cout << endl;

    vector <vector <double> > unit_cell_xyz = get_simplecubic(3, 3, 3, 1);
    cout << "NEIGHBOUR LIST FUNCTION TESTING PART" << endl;
    cout << "Atomic position of unit cell, simple cubic" << endl;
    print_2dvector(unit_cell_xyz);
    double a = 3;
    vector <vector <double> > unit_cell_vec = {{a, 0, 0}, {0, a, 0}, {0, 0, a}};
    cout << "Unit cell vectors" << endl;
    print_2dvector(unit_cell_vec);
    // Call neighbout list generation with cutoff 1.1 Angstrom (Nearest separation of my simple cubic case is 1 angstrom)
    write_neighbour_list(neighbour_list(unit_cell_vec, unit_cell_xyz, 1.1));
    cout << "Output neighbout list in csv" << endl;
    return 0;
}