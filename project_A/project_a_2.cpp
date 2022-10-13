#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
// Project A part 1 functions are reused
// simple cubic, bcc, fcc, diamond structure generation
#include "project_a_1.h"
#include "project_a_2.h"
using namespace std;

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

    // vector <vector <double> > unit_cell_xyz = get_simplecubic(3, 3, 3, 1);
    vector <vector <double> > unit_cell_xyz = get_diamond(3, 3, 3, 1);
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