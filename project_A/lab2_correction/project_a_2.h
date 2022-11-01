#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;

// Custom struct for storing neighbour records
struct neighbour_record {
    int label1;
    int label2;
    double squared_distance;
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
        sum += u[i]*v[i];
    return sum;
}

// Volume (3-D)
double volume(vector <vector <double> > lattice) {
    return dot_prod(lattice[0], cross_prod(lattice[1], lattice[2]));
}

// Reciprocal vector calculation
// input 3 lattice vectors, output 3 reciprocal lattice vectors
vector <vector <double> > reciprocal_vec(vector <vector <double> > a) {
    vector <double> crossprod_temp = cross_prod(a[1], a[2]);
    double volume = dot_prod(a[0], crossprod_temp);
    
    vector <double> b1(3);
    for (int i = 0; i < crossprod_temp.size(); i++)
        b1[i] = (crossprod_temp[i] / volume);

    crossprod_temp = cross_prod(a[2], a[0]);
    vector <double> b2(3);
    for (int i = 0; i < crossprod_temp.size(); i++)
        b2[i] = (crossprod_temp[i] / volume);

    crossprod_temp = cross_prod(a[0], a[1]);
    vector <double> b3(3);
    for (int i = 0; i < crossprod_temp.size(); i++)
        b3[i] = (crossprod_temp[i] / volume);

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
        print_1dvector(v);
    return;
}

// PBC on fractional coords
// input a 1-D vector, output a 1-D vector with coefficients in [-0.5, 0.5]
// vector <double> pbc(vector <double> d) {
//     double temp;
//     vector <double> res(3);
//     for (int i = 0; i < d.size(); i++)
//     {
//         temp = fmod(d[i], 1.0);
//         if (temp > 0.5)
//         {
//             temp -= 1;
//         }
//         else if (temp < -0.5)
//         {
//             temp += 1;
//         }
//         res[i] = temp;
//     }
//     return res;
// }
double pbc(double d) {
    double temp;
    vector <double> res(3);
    temp = fmod(d, 1.0);
    if (temp > 0.5)
        temp -= 1;
    else if (temp < -0.5)
        temp += 1;
    return temp;
}

// True distance
// input: distance vector (1-D), real lattice vectors, reciprocal lattice vectors, print option
// output: PBC reciprocal distance vector and PBC real space distance vector
pair<vector <double>, vector <double>>true_distance_vec(vector <double> d, vector <vector <double> > real_lattice, vector <vector <double> > recip_lattice, bool print=false) {
    vector <double> recip_d(3), true_d(3);
    if (print)
        cout << "Reciprocal distance vector" << endl;
    // for (auto recip_vec: recip_lattice)
    // {
    //     // Step 1: Compute reciprocal distance vector
    //     // Each reciprocal lattice vector dot distance vector (1D)
    //     // Yields reciprocal projection of distance vector
    //     recip_d.push_back(dot_prod(recip_vec, d));
    //     if (print)
    //         cout << recip_d.back() << " ";
    // }
    for (int i = 0; i < recip_lattice.size(); i++)
        recip_d[i] = pbc(dot_prod(recip_lattice[i], d));
    // Step 2: Apply PBC on reciprocal distance vector
    // i.e. PBC in reciprocal space
    // recip_d = pbc(recip_d);
    if (print) {
        cout << endl;
        cout << "PBC reciprocal distance vector" << endl;
        print_1dvector(recip_d);
        cout << "PBC real distance vector" << endl;
    }
    // for (auto real_vec: real_lattice)
    // {
    //     // Step 3: Convert the PBC reciprocal distance vector to real space
    //     // Each real lattice vector dot PBC reciprocal distance vector
    //     // Yields PBC real space distance vector
    //     true_d.push_back(dot_prod(real_vec, recip_d));
    //     if (print)
    //         cout << true_d.back() << " ";
    // }
    for (int i = 0; i < real_lattice.size(); i++)
    {
        true_d[i] = dot_prod(real_lattice[i], recip_d);
        if (print)
            cout << true_d.back() << " ";
    }

    if (print)
        cout << endl;
    // first stores reciprocal distance, second stores true distance
    return make_pair(recip_d, true_d);
}

// Neighbout list function
// input: real lattice vectors, xyz coords of all atoms in the structure, distance cutoff
// output: neighbour list (Array of neighbour record struct)
vector <neighbour_record> neighbour_list(vector <vector <double> > lattice, vector <vector <double> > xyz, double dist_cutoff) {
    vector <neighbour_record> nlist;
    nlist.reserve(xyz.size() * xyz.size());
    vector <double> dist_vec_temp, true_dist_vec_temp;
    double true_dist_temp;
    int idx = 0;
    // Iterate all possible pair of atoms (without self-self counting and double counting)
    for (int i = 0; i < xyz.size(); i++)
    {
        for (int j = i + 1; j < xyz.size(); j++)
        {
            // Step 1: Compute distance vector of chosen pair of atoms
            dist_vec_temp = {xyz[j][0] - xyz[i][0], xyz[j][1] - xyz[i][1], xyz[j][2] - xyz[i][2]};
            // Step 2: Get true distance vector in real space (PBC applied)
            true_dist_vec_temp = true_distance_vec(dist_vec_temp, lattice, reciprocal_vec(lattice)).second;
            // Step 3: Self dot product of true distance vector to get distance squared
            true_dist_temp = dot_prod(true_dist_vec_temp, true_dist_vec_temp);
            // Step 4: if distance squared is smaller than distance cutoff squared, add the pair to neighbour list
            if (true_dist_temp < (dist_cutoff * dist_cutoff))
            {
                // cout << i << " " << j << " " << true_dist_temp << " " << (dist_cutoff * dist_cutoff) << endl;
                nlist.push_back({i, j, true_dist_temp});
                // ### 
                // Correcion here, j is also neighbour of i 
                // ###
                nlist.push_back({j, i, true_dist_temp});
            }       
        }
    }
    return nlist;
}

// Write the neighbour list into a csv file
void write_neighbour_list(vector <neighbour_record> nlist, string filename="neighbour_list.csv") {
    ofstream file;
    file.open(filename);
    file << "label1,label2,squared_distance\n";
    for (auto record: nlist) {
        file << record.label1 << "," << record.label2 << "," << record.squared_distance << "\n";
    }
    file.close();
    return;
}