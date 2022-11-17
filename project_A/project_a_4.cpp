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
vector <vector <double> > add_perturbation(vector <vector <double> > xyz, double perturb_strength=0.01) {
    srand (time(NULL));
    // srand (1);
    for (int i = 0; i < xyz.size(); i++) {
        for (int j = 0; j < 3; j++) {
            xyz[i][j] += ((rand() % 100000) / 50000. - 1) * 0.01;
        }
    }
    return xyz;
}

// Q is the generalized coordinates, e.g. dx or dy or dz
// double LJ_slope(double r2, double Q, double epsilon, double sigma) {
double LJ_slope_prefactor(double r2, double epsilon, double sigma) {
    double sigma2 = sigma * sigma;
    double s_r_6 = (sigma2 / r2) * (sigma2 / r2) * (sigma2 / r2);
    return 4 * epsilon * (- 12 * s_r_6 * s_r_6 / r2 + 6 * s_r_6 / r2);
}

double CoulombBuckingham_slope_prefactor(double r2, double q1, double q2, double A, double B, double C) {
    double r = sqrt(r2);
    return -A * B * exp(-B * r) / r + 6 * C / (r2 * r2 * r2 * r2) + q1 * q2 / (4 * pi * e_0 * r2 * r);
}

double total_LJ(vector <neighbour_record> nlist, double epsilon, double sigma) {
    double en_tot = 0;
    for (auto record: nlist) {
        en_tot += LJ(record.squared_distance, epsilon, sigma);
    }
    return en_tot;
}

double total_CoulombBuckingham(vector <neighbour_record> nlist, int length) {
    double en_tot = 0;
    for (auto record: nlist) {
        if ((record.label1 < length) && (record.label2 < length)) {
            en_tot += CoulombBuckingham(record.squared_distance, q_Mg, q_Mg, A_MgMg, B_MgMg, C_MgMg);
        }
        // Both ions are O
        else if ((record.label1 >= length) && (record.label2 >= length)) {
            en_tot += CoulombBuckingham(record.squared_distance, q_O, q_O, A_OO, B_OO, C_OO);
        }
        // One Mg One O
        else {
            en_tot += CoulombBuckingham(record.squared_distance, q_Mg, q_O, A_MgO, B_MgO, C_MgO);
        }
    }
    return en_tot;
}

vector <vector <double> > Fprime_CoulombBuckingham(vector <neighbour_record> nlist, int length, int Mg_length) {
    vector <vector <double> > Fprime(length, vector<double>(3, 0));
    double slope_pref;
    for (auto record: nlist) {
        // Computing the gradient vector elements for each pair
        // Summing all the neighbour effects
        // The first atom (vector end) get the gradient
        if ((record.label1 < Mg_length) && (record.label2 < Mg_length)) {
            slope_pref = CoulombBuckingham_slope_prefactor(record.squared_distance, q_Mg, q_Mg, A_MgMg, B_MgMg, C_MgMg);
        }
        // Both ions are O
        else if ((record.label1 >= Mg_length) && (record.label2 >= Mg_length)) {
            slope_pref = CoulombBuckingham_slope_prefactor(record.squared_distance, q_O, q_O, A_OO, B_OO, C_OO);
        }
        // One Mg One O
        else {
            slope_pref = CoulombBuckingham_slope_prefactor(record.squared_distance, q_Mg, q_O, A_MgO, B_MgO, C_MgO);
        }
        
        Fprime[record.label1][0] -= slope_pref * record.dx;
        Fprime[record.label1][1] -= slope_pref * record.dy;
        Fprime[record.label1][2] -= slope_pref * record.dz;
        // The second atom (vector head) get the opposite gradient
        Fprime[record.label2][0] += slope_pref * record.dx;
        Fprime[record.label2][1] += slope_pref * record.dy;
        Fprime[record.label2][2] += slope_pref * record.dz;
    }
    return Fprime;
}

double secant_method(vector <double> Fprime, vector <double> h, vector <double> Fprime_step, double rate=0.1, double tol=1e-5) {
    double alpha;
    vector <double> diff(3);
    diff[0] = Fprime_step[0] - Fprime[0];
    diff[1] = Fprime_step[1] - Fprime[1];
    diff[2] = Fprime_step[2] - Fprime[2];
    if (dot_prod(diff, diff) < tol) {
        return 0.;
    } 
    alpha = -rate * dot_prod(Fprime, h) / dot_prod(diff, h);
    return alpha;
}

void SD_LJ(vector <vector <double> > &xyz, vector <vector <double> > lat_vec, double rate=0.1, double grad_tol=1e-5, double nn_cutoff=numeric_limits<int>::max()) {
    vector <vector <double> > Fprime(xyz.size(), vector<double>(3, 0));
    vector <neighbour_record> nlist = neighbour_list(lat_vec, xyz, nn_cutoff);
    // double Fprime_xtmp, Fprime_ytmp, Fprime_ztmp;
    vector <double> Fprime_tmp(3, 0);
    double slope_pref;
    // print_2dvector(xyz);
    // Compute full gradient vectors for all atoms
    for (auto record: nlist) {
        // Computing the gradient vector elements for each pair
        // Summing all the neighbour effects
        // The first atom (vector end) get the gradient
        slope_pref = LJ_slope_prefactor(record.squared_distance, epsilon_ArAr, sigma_ArAr);
        Fprime[record.label1][0] -= slope_pref * record.dx;
        Fprime[record.label1][1] -= slope_pref * record.dy;
        Fprime[record.label1][2] -= slope_pref * record.dz;
        // The second atom (vector head) get the opposite gradient
        Fprime[record.label2][0] += slope_pref * record.dx;
        Fprime[record.label2][1] += slope_pref * record.dy;
        Fprime[record.label2][2] += slope_pref * record.dz;
    }

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
        
        slope_pref = LJ_slope_prefactor(step_r2_tmp, epsilon_ArAr, sigma_ArAr);

        Fprime_step[record.label1][0] -= slope_pref * dQ_tmp[0];
        Fprime_step[record.label1][1] -= slope_pref * dQ_tmp[1];
        Fprime_step[record.label1][2] -= slope_pref * dQ_tmp[2];

        Fprime_step[record.label2][0] += slope_pref * dQ_tmp[0];
        Fprime_step[record.label2][1] += slope_pref * dQ_tmp[1];
        Fprime_step[record.label2][2] += slope_pref * dQ_tmp[2];
    }

    double alpha;
    vector <double> diff(3), h(3);
    for (int n = 0; n < xyz.size(); n++) {
        h[0] = -Fprime[n][0];
        h[1] = -Fprime[n][1];
        h[2] = -Fprime[n][2];
        alpha = secant_method(Fprime[n], h, Fprime_step[n], rate, grad_tol);
        // cout << alpha << endl;
        xyz[n][0] += alpha * -Fprime[n][0];
        xyz[n][1] += alpha * -Fprime[n][1];
        xyz[n][2] += alpha * -Fprime[n][2];
    }
}

void CG_ColumbBuckingham(vector <vector <double> > &xyz, vector <vector <double> > lat_vec, vector <vector <double> > &h, int iter, double rate=0.1, double grad_tol=1e-5, double nn_cutoff=numeric_limits<int>::max()) {
    int length = xyz.size();
    vector <vector <double> > Fprime(length, vector<double>(3, 0));
    vector <neighbour_record> nlist = neighbour_list(lat_vec, xyz, nn_cutoff);
    // double Fprime_xtmp, Fprime_ytmp, Fprime_ztmp;
    vector <double> Fprime_tmp(3, 0);
    double slope_pref;
    // print_2dvector(xyz);
    // Compute full gradient vectors for all atoms
    Fprime = Fprime_CoulombBuckingham(nlist, length, length / 2);

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
        
        if ((record.label1 < length) && (record.label2 < length)) {
            slope_pref = CoulombBuckingham_slope_prefactor(step_r2_tmp, q_Mg, q_Mg, A_MgMg, B_MgMg, C_MgMg);
        }
        // Both ions are O
        else if ((record.label1 >= length) && (record.label2 >= length)) {
            slope_pref = CoulombBuckingham_slope_prefactor(step_r2_tmp, q_O, q_O, A_OO, B_OO, C_OO);
        }
        // One Mg One O
        else {
            slope_pref = CoulombBuckingham_slope_prefactor(step_r2_tmp, q_Mg, q_O, A_MgO, B_MgO, C_MgO);
        }

        Fprime_step[record.label1][0] -= slope_pref * dQ_tmp[0];
        Fprime_step[record.label1][1] -= slope_pref * dQ_tmp[1];
        Fprime_step[record.label1][2] -= slope_pref * dQ_tmp[2];

        Fprime_step[record.label2][0] += slope_pref * dQ_tmp[0];
        Fprime_step[record.label2][1] += slope_pref * dQ_tmp[1];
        Fprime_step[record.label2][2] += slope_pref * dQ_tmp[2];
    }

    double alpha;
    vector <double> diff(3), h_new(3);
    for (int n = 0; n < xyz.size(); n++) {
        alpha = secant_method(Fprime[n], h[n], Fprime_step[n], rate, grad_tol);
        // cout << alpha << endl;
        xyz[n][0] += alpha * h[n][0];
        xyz[n][1] += alpha * h[n][1];
        xyz[n][2] += alpha * h[n][2];
    }

    nlist = neighbour_list(lat_vec, xyz, nn_cutoff);
    vector <vector <double> > Fprime_new = Fprime_CoulombBuckingham(nlist, length, length / 2);
    for (int n = 0; n < xyz.size(); n++) {
        double gamma = dot_prod(Fprime_new[n], Fprime_new[n]) / dot_prod(Fprime[n], Fprime[n]);
        h[n][0] += -Fprime[n][0] + gamma * h[n][0];
        h[n][1] += -Fprime[n][1] + gamma * h[n][1];
        h[n][2] += -Fprime[n][2] + gamma * h[n][2];
    }
}

int main() {
    vector <vector <double> > Ar_cell_xyz, Ar_cell_vec;
    int dim = 2;
    int n_cells = dim;
    double en, en_tmp, diff;
    double en_tol = 1e-4;
    int n_iter = 20000;
    Ar_cell_xyz = get_fcc(n_cells, n_cells, n_cells, a_ArAr);
    Ar_cell_vec = {{n_cells * a_ArAr, 0, 0}, {0, n_cells * a_ArAr, 0}, {0, 0, n_cells * a_ArAr}};
    Ar_cell_xyz = add_perturbation(Ar_cell_xyz, 0.01);

    // cout << "Perturbed Structure Coordinates" << endl;
    // print_2dvector(Ar_cell_xyz);
    // double nn_cutoff = 3.1*a_ArAr;
    double nn_cutoff = numeric_limits<int>::max();
    en_tmp = total_LJ(neighbour_list(Ar_cell_vec, Ar_cell_xyz, nn_cutoff), epsilon_ArAr, sigma_ArAr);
    cout << "Initial energy " << en_tmp << endl;
    for (int i = 0; i < n_iter; i++)
    {
        SD_LJ(Ar_cell_xyz, Ar_cell_vec, 0.001, 1e-5, nn_cutoff);
        if ((i+1) % 500 == 0) {
            en = total_LJ(neighbour_list(Ar_cell_vec, Ar_cell_xyz, nn_cutoff), epsilon_ArAr, sigma_ArAr);
            diff = en - en_tmp;
            cout << "Step " << i+1 << "; energy: " << en << "; abs difference: " << fabs(diff) << "; convergence tol: " << en_tol << endl;
            if (fabs(diff) < en_tol) {
                cout << "convergence reached, terminate SD" << endl;
                break;
            }
            en_tmp = en;
        }
    }
    // cout << "Relaxed Structure Coordinates" << endl;
    // print_2dvector(Ar_cell_xyz);

    vector <vector <double> > MgO_cell_vec, MgO_cell_xyz, O_cell_xyz;
    MgO_cell_vec = {{n_cells * a_MgO, 0, 0}, {0, n_cells * a_MgO, 0}, {0, 0, n_cells * a_MgO}};
    MgO_cell_xyz = get_fcc(n_cells, n_cells, n_cells, a_MgO);
    int Mg_length = MgO_cell_xyz.size(); 
    O_cell_xyz = get_fcc(n_cells, n_cells, n_cells, a_MgO, 0.5);
    MgO_cell_xyz.insert(end(MgO_cell_xyz), begin(O_cell_xyz), end(O_cell_xyz));

    MgO_cell_xyz = add_perturbation(MgO_cell_xyz, 0.01);

    vector <neighbour_record> MgO_nlist = neighbour_list(MgO_cell_vec, MgO_cell_xyz, numeric_limits<int>::max());
    vector <vector <double> > Fprime_0 = Fprime_CoulombBuckingham(MgO_nlist, MgO_cell_xyz.size(), Mg_length);
    vector <vector <double> > h(MgO_cell_xyz.size(), vector<double>(3, 0));
    for (int i = 0; i < MgO_cell_xyz.size(); i++) {
        h[i][0] = -Fprime_0[i][0];
        h[i][1] = -Fprime_0[i][1];
        h[i][2] = -Fprime_0[i][2];
    }

    en_tmp = total_CoulombBuckingham(neighbour_list(MgO_cell_vec, MgO_cell_xyz, nn_cutoff), Mg_length);
    cout << "Initial energy " << en_tmp << endl;
    for (int i = 0; i < n_iter; i++)
    {
        CG_ColumbBuckingham(MgO_cell_xyz, MgO_cell_vec, h, i, 0.1, 1e-5, nn_cutoff);
        if ((i+1) % 500 == 0) {
            en = total_CoulombBuckingham(neighbour_list(MgO_cell_vec, MgO_cell_xyz, nn_cutoff), Mg_length);
            diff = en - en_tmp;
            cout << "Step " << i+1 << "; energy: " << en << "; abs difference: " << fabs(diff) << "; convergence tol: " << en_tol << endl;
            if (fabs(diff) < en_tol) {
                cout << "convergence reached, terminate SD" << endl;
                break;
            }
            en_tmp = en;
        }
    }
    return 0;
}