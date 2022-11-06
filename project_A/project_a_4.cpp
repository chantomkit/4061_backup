#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <limits>
#include <time.h>
#include <stdlib.h>
// simple cubic, bcc, fcc, diamond structure generation
#include "project_a_1.h"
// for generating neighbour list for 2 body potentials
#include "project_a_2.h"
// Constants and potentials
#include "project_a_3.h"
using namespace std;

// add random perturb (-0.1 ~ 0.1)
vector <vector <double> > add_perturbation(vector <vector <double> > xyz) {
    srand (time(NULL));
    for (int i = 0; i < xyz.size(); i++) {
        for (int j = 0; j < 3; j++) {
            xyz[i][j] += ((rand() % 100000) / 50000. - 1) * 0.001;
        }
    }
    return xyz;
}

// Q is the generalized coordinates, e.g. dx or dy or dz
double LJ_slope(double r2, double Q, double epsilon, double sigma) {
    double r = sqrt(r2);
    double sq_sigma = sigma * sigma;
    double s_r_6 = (sq_sigma / r2) * (sq_sigma / r2) * (sq_sigma / r2);
    return 4 * epsilon * Q * (- 12 * s_r_6 * s_r_6 / r2 + 6 * s_r_6 / r2);
}

vector <vector <double> > SD(vector <vector <double> > xyz, vector <vector <double> > lat_vec, double rate=0.01) {
    vector <vector <double> > Fprime(xyz.size(), vector<double>(3, 0));
    vector <neighbour_record> nlist = neighbour_list(lat_vec, xyz, numeric_limits<int>::max());
    // double Fprime_xtmp, Fprime_ytmp, Fprime_ztmp;
    vector <double> Fprime_tmp(3, 0);

    // print_2dvector(xyz);
    // Compute full gradient vectors for all atoms
    for (auto record: nlist) {
        // Computing the gradient vector elements for each pair
        // Summing all the neighbour effects
        Fprime_tmp[0] = LJ_slope(record.squared_distance, record.dx, epsilon_ArAr, sigma_ArAr);
        Fprime_tmp[1] = LJ_slope(record.squared_distance, record.dy, epsilon_ArAr, sigma_ArAr);
        Fprime_tmp[2] = LJ_slope(record.squared_distance, record.dz, epsilon_ArAr, sigma_ArAr);
        // The first atom (vector end) get the gradient
        Fprime[record.label1][0] -= Fprime_tmp[0];
        Fprime[record.label1][1] -= Fprime_tmp[1];
        Fprime[record.label1][2] -= Fprime_tmp[2];
        // The second atom (vector head) get the opposite gradient
        Fprime[record.label2][0] += Fprime_tmp[0];
        Fprime[record.label2][1] += Fprime_tmp[1];
        Fprime[record.label2][2] += Fprime_tmp[2];
    }
    // print_2dvector(Fprime);

    vector <vector <double> > Fprime_step(xyz.size(), vector<double>(3, 0));
    vector <double> step1_tmp(3), step2_tmp(3), dQ_tmp(3);
    double step_r2_tmp;
    // Compute gradient vectors of all atoms after taking small step
    for (auto record: nlist) {
        step1_tmp[0] = xyz[record.label1][0] + rate * -Fprime[record.label1][0];
        step1_tmp[1] = xyz[record.label1][1] + rate * -Fprime[record.label1][1];
        step1_tmp[2] = xyz[record.label1][2] + rate * -Fprime[record.label1][2];

        step2_tmp[0] = xyz[record.label2][0] + rate * -Fprime[record.label2][0];
        step2_tmp[1] = xyz[record.label2][1] + rate * -Fprime[record.label2][1];
        step2_tmp[2] = xyz[record.label2][2] + rate * -Fprime[record.label2][2];

        dQ_tmp[0] = step2_tmp[0] - step1_tmp[0];
        dQ_tmp[1] = step2_tmp[1] - step1_tmp[1];
        dQ_tmp[2] = step2_tmp[2] - step1_tmp[2];
        step_r2_tmp = dot_prod(dQ_tmp, dQ_tmp);
        
        Fprime_tmp[0] = LJ_slope(step_r2_tmp, dQ_tmp[0], epsilon_ArAr, sigma_ArAr);
        Fprime_tmp[1] = LJ_slope(step_r2_tmp, dQ_tmp[1], epsilon_ArAr, sigma_ArAr);
        Fprime_tmp[2] = LJ_slope(step_r2_tmp, dQ_tmp[2], epsilon_ArAr, sigma_ArAr);

        Fprime_step[record.label1][0] -= Fprime_tmp[0];
        Fprime_step[record.label1][1] -= Fprime_tmp[1];
        Fprime_step[record.label1][2] -= Fprime_tmp[2];

        Fprime_step[record.label2][0] += Fprime_tmp[0];
        Fprime_step[record.label2][1] += Fprime_tmp[1];
        Fprime_step[record.label2][2] += Fprime_tmp[2];

        // step_tmp[0] = xyz[record.label2][0] + rate * -Fprime[record.label2][0];
        // step_tmp[1] = xyz[record.label2][1] + rate * -Fprime[record.label2][1];
        // step_tmp[2] = xyz[record.label2][2] + rate * -Fprime[record.label2][2];
        // dQ_tmp[0] = xyz[record.label1][0] - step_tmp[0];
        // dQ_tmp[1] = xyz[record.label1][1] - step_tmp[1];
        // dQ_tmp[2] = xyz[record.label1][2] - step_tmp[2];
        
        // Fprime_step[record.label2][0] += LJ_slope(dot_prod(dQ_tmp, dQ_tmp), dQ_tmp[0], epsilon_ArAr, sigma_ArAr);
        // Fprime_step[record.label2][1] += LJ_slope(dot_prod(dQ_tmp, dQ_tmp), dQ_tmp[1], epsilon_ArAr, sigma_ArAr);
        // Fprime_step[record.label2][2] += LJ_slope(dot_prod(dQ_tmp, dQ_tmp), dQ_tmp[2], epsilon_ArAr, sigma_ArAr);
    }
    // cout << endl;
    // print_2dvector(Fprime_step);

    double alpha;
    vector <double> diff(3), h(3);
    for (int n = 0; n < xyz.size(); n++) {
        diff[0] = Fprime_step[n][0] - Fprime[n][0];
        diff[1] = Fprime_step[n][1] - Fprime[n][1];
        diff[2] = Fprime_step[n][2] - Fprime[n][2];
        h[0] = -Fprime[n][0];
        h[1] = -Fprime[n][1];
        h[2] = -Fprime[n][2];
        alpha = -rate * dot_prod(Fprime[n], h) / dot_prod(diff, h);
        if (alpha > 1e-6) {
            xyz[n][0] += alpha * -Fprime[n][0];
            xyz[n][1] += alpha * -Fprime[n][1];
            xyz[n][2] += alpha * -Fprime[n][2];
        } 
    }
    // cout << "\n" << endl;
    // print_2dvector(xyz);
    return xyz;
}

int main() {
    vector <vector <double> > Ar_cell_xyz, Ar_cell_vec;
    int dim = 1;
    int n_cells = dim;
    Ar_cell_xyz = get_fcc(n_cells, n_cells, n_cells, a_ArAr);
    Ar_cell_vec = {{n_cells * a_ArAr, 0, 0}, {0, n_cells * a_ArAr, 0}, {0, 0, n_cells * a_ArAr}};
    write_xyz(Ar_cell_xyz);
    write_neighbour_list(neighbour_list(Ar_cell_vec, Ar_cell_xyz, 10));
    
    Ar_cell_xyz = add_perturbation(Ar_cell_xyz);
    print_2dvector(Ar_cell_xyz);
    cout << endl;
    for (int i = 0; i < 200; i++)
    {
        Ar_cell_xyz = SD(Ar_cell_xyz, Ar_cell_vec);
        // break;
    }
    cout << endl;  
    print_2dvector(Ar_cell_xyz);
    return 0;
}