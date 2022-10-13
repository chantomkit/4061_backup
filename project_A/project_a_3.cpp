#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
// Project A part 1 functions are reused
// simple cubic, bcc, fcc, diamond structure generation
#include "project_a_1.h"
#include "project_a_2.h"
using namespace std;

// Argon LJ params: https://www.researchgate.net/figure/Lennard-Jones-LJ-potential-parameters-of-different-materials-considered-in-thepresent_tbl2_319412425
// Argon structure and lattice constant: https://aip.scitation.org/doi/10.1063/1.1726009?cookieSet=1 
// Argon cohesive energy: https://arxiv.org/pdf/2012.05413.pdf 
// https://www.knowledgedoor.com/2/elements_handbook/cohesive_energy.html
//       | epsilon (eV) | sigma (nm) | a (Angstrom) | structure | cohesive energy (eV) |
// Ar-Ar | 0.34         | 0.0104     | 5.311        | FCC       | 0.08                 |
const double pi = 3.141592654;
const double e_0 = 55.26349406 / 1e4; // e^2 eV^-1 Angstrom^-1
const double eV = 1.602176634e-19;

const double a_ArAr = 5.311; // Angstrom
const double epsilon_ArAr = 0.0104; // eV
const double sigma_ArAr = 0.34 * 10; // Angstrom

const double En_Ar = -0.08; // eV

// MgO Buckingham params: https://www.researchgate.net/publication/248076144_Effect_of_inversion_on_thermoelastic_and_thermal_transport_properties_of_MgAl_2_O_4_spinel_by_atomistic_simulation
// MgO structural properties: https://materialsproject.org/materials/mp-1265/
// MgO Cohesive energy: https://iopscience.iop.org/article/10.1088/1742-6596/377/1/012067/pdf 
//      | a (Angstrom) | structure          | cohesive energy (eV) |
// Mg-O | 4.26         | FCC (2 atom basis) | 10.98                |
const double a_MgO = 4.26; // Angstrom
const double A_MgO = 1279.69; // eV
const double B_MgO = 1. / 0.2997; // Angstrom^-1
const double C_MgO = 0; // eV Angstrom^6
const double q_Mg = 2; // e
const double q_O = -2; // e

const double En_MgO = -10.98; // eV

double LJ(double r2, double epsilon, double sigma) {
    double sq_sigma = sigma * sigma;
    double s_r_6 = (sq_sigma / r2) * (sq_sigma / r2) * (sq_sigma / r2);
    return 4 * epsilon * (s_r_6 * s_r_6 - s_r_6);
}

double CoulombBuckingham(double r2, double q1, double q2, double A, double B, double C) {
    double r = sqrt(r2);
    return A * exp(-B * r) - (C / (r2 * r2 * r2)) + (q1 * q2) / (4 * pi * e_0 * r);
}

int main() {
    int n_cells = 3;
    vector <vector <double> > Ar_cell_xyz = get_fcc(n_cells, n_cells, n_cells, a_ArAr);
    cout << n_cells << " by " << n_cells << " Cell" << endl << "N atoms: " << Ar_cell_xyz.size() << endl;
    vector <vector <double> > Ar_cell_vec = {{n_cells * a_ArAr, 0, 0}, {0, n_cells * a_ArAr, 0}, {0, 0, n_cells * a_ArAr}};
    vector <neighbour_record> nlist = neighbour_list(Ar_cell_vec, Ar_cell_xyz, 5.1 * a_ArAr);

    double tot_en = 0;
    cout << "LJ potential" << endl;
    for (auto record: nlist) {
        tot_en += LJ(record.squared_distance, epsilon_ArAr, sigma_ArAr);
        // cout << record.label1 << " " << record.label2 << " " << record.squared_distance << " " << tot_en << endl;
    }
    cout << "Computed Total Cohesive Energy of Ar (eV): " << tot_en << endl << "Computed Energy per Ar atom (eV): " << tot_en / Ar_cell_xyz.size() << endl;
    cout << "Reference value of Cohesive Energy of Ar per atom (eV): " << En_Ar << endl << endl;

    vector <vector <double> > MgO_cell_vec = {{n_cells * a_MgO, 0, 0}, {0, n_cells * a_MgO, 0}, {0, 0, n_cells * a_MgO}};
    // The name is MgO cell and MgO charges but firstly initialize Mg atoms
    vector <vector <double> > MgO_cell_xyz = get_fcc(n_cells, n_cells, n_cells, a_MgO);
    vector <double> MgO_charges(MgO_cell_xyz.size(), q_Mg);
    // Then initialize O atoms
    vector <vector <double> > O_cell_xyz = get_fcc(n_cells, n_cells, n_cells, a_MgO, a_MgO/2);
    vector <double> O_charges(O_cell_xyz.size(), q_O);
    // Insert O atoms to the end of Mg atoms to make MgO cell
    MgO_cell_xyz.insert(end(MgO_cell_xyz), begin(O_cell_xyz), end(O_cell_xyz));
    MgO_charges.insert(end(MgO_charges), begin(O_charges), end(O_charges));
    vector <neighbour_record> MgO_nlist = neighbour_list(MgO_cell_vec, MgO_cell_xyz, 5.1 * a_ArAr);
    double MgO_tot_en = 0;
    double q1, q2;
    cout << n_cells << " by " << n_cells << " Cell" << endl << "N atoms: " << MgO_cell_xyz.size() << endl;
    cout << "Coulomb-Buckingham Potential" << endl;
    for (auto record: MgO_nlist) {
        q1 = MgO_charges[record.label1];
        q2 = MgO_charges[record.label2];
        MgO_tot_en += CoulombBuckingham(record.squared_distance, q1, q2, A_MgO, B_MgO, C_MgO);
        // cout << record.label1 << " " << record.label2 << " " << record.squared_distance << " " << tot_en << endl;
    }
    cout << "Computed Total Cohesive Energy of MgO (eV): " << MgO_tot_en << endl << "Computed Energy per MgO atom (eV): " << MgO_tot_en / MgO_cell_xyz.size() << endl;
    cout << "Reference value of Cohesive Energy of MgO per atom (eV): " << En_MgO << endl;
}