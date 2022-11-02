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
            xyz[i][j] += ((rand() % 100000) / 50000. - 1) * 0.1;
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

vector <vector <double> > SD(vector <vector <double> > xyz, vector <vector <double> > lat_vec, double rate=1e-6) {
    vector <vector <double> > Fprime(xyz.size(), vector<double>(3, 0));
    vector <neighbour_record> nlist = neighbour_list(lat_vec, xyz, numeric_limits<int>::max());
    double Fprime_xtmp, Fprime_ytmp, Fprime_ztmp;

    // print_2dvector(xyz);
    // Compute full gradient vectors for all atoms
    for (auto record: nlist) {
        // Computing the gradient vector elements for each pair
        // Summing all the neighbour effects
        Fprime_xtmp = LJ_slope(record.squared_distance, record.dx, epsilon_ArAr, sigma_ArAr);
        Fprime_ytmp = LJ_slope(record.squared_distance, record.dy, epsilon_ArAr, sigma_ArAr);
        Fprime_ztmp = LJ_slope(record.squared_distance, record.dz, epsilon_ArAr, sigma_ArAr);
        // The first atom (vector end) get the gradient
        Fprime[record.label1][0] += Fprime_xtmp;
        Fprime[record.label1][1] += Fprime_ytmp;
        Fprime[record.label1][2] += Fprime_ztmp;
        // The second atom (vector head) get the opposite gradient
        Fprime[record.label2][0] -= Fprime_xtmp;
        Fprime[record.label2][1] -= Fprime_ytmp;
        Fprime[record.label2][2] -= Fprime_ztmp;
    }

    vector <vector <double> > Fprime_step(xyz.size(), vector<double>(3, 0));
    vector <double> step_tmp(3);
    double step_r2;
    // Compute full gradient vectors for all atoms
    for (auto record: nlist) {
        // Computing the gradient vector elements for each pair
        // Summing all the neighbour effects
        // The first atom (vector end) get the gradient
        step_tmp[0] = record.dx + rate * -Fprime[record.label1][0];
        step_tmp[1] = record.dy + rate * -Fprime[record.label1][1];
        step_tmp[2] = record.dz + rate * -Fprime[record.label1][2];
        step_r2 = dot_prod(step_tmp, step_tmp);
        
        Fprime_step[record.label1][0] += LJ_slope(step_r2, step_tmp[0], epsilon_ArAr, sigma_ArAr);;
        Fprime_step[record.label1][1] += LJ_slope(step_r2, step_tmp[1], epsilon_ArAr, sigma_ArAr);;
        Fprime_step[record.label1][2] += LJ_slope(step_r2, step_tmp[2], epsilon_ArAr, sigma_ArAr);;

        // The second atom (vector head) get the opposite gradient
        step_tmp[0] = record.dx + rate * -Fprime[record.label2][0];
        step_tmp[1] = record.dy + rate * -Fprime[record.label2][1];
        step_tmp[2] = record.dz + rate * -Fprime[record.label2][2];
        step_r2 = dot_prod(step_tmp, step_tmp);
        
        Fprime_step[record.label2][0] += LJ_slope(step_r2, -step_tmp[0], epsilon_ArAr, sigma_ArAr);;
        Fprime_step[record.label2][1] += LJ_slope(step_r2, -step_tmp[1], epsilon_ArAr, sigma_ArAr);;
        Fprime_step[record.label2][2] += LJ_slope(step_r2, -step_tmp[2], epsilon_ArAr, sigma_ArAr);;
    }

    double alpha;
    vector <double> diff(3);
    for (int n = 0; n < xyz.size(); n++) {
        diff[0] = Fprime_step[n][0] - Fprime[n][0];
        diff[1] = Fprime_step[n][1] - Fprime[n][1];
        diff[2] = Fprime_step[n][2] - Fprime[n][2];
        alpha = -rate * dot_prod(Fprime[n], Fprime[n]) / dot_prod(diff, Fprime[n]);
        
        xyz[n][0] += alpha * -Fprime[n][0];
        xyz[n][1] += alpha * -Fprime[n][1];
        xyz[n][2] += alpha * -Fprime[n][2];
    }
    // cout << "\n" << endl;
    // print_2dvector(xyz);
    return xyz;
}

int main() {
    vector <vector <double> > Ar_cell_xyz, Ar_cell_vec;
    int dim = 1;
    int n_cells = dim*2;
    Ar_cell_xyz = get_fcc(n_cells, n_cells, n_cells, a_ArAr);
    Ar_cell_vec = {{n_cells * a_ArAr, 0, 0}, {0, n_cells * a_ArAr, 0}, {0, 0, n_cells * a_ArAr}};
    
    Ar_cell_xyz = add_perturbation(Ar_cell_xyz);
    print_2dvector(Ar_cell_xyz);
    for (int i = 0; i < 10; i++)
        Ar_cell_xyz = SD(Ar_cell_xyz, Ar_cell_vec);
    cout << "\n" << endl;
    print_2dvector(Ar_cell_xyz);
    return 0;
}