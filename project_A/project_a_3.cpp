#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
// Project A part 1 functions are reused
// simple cubic, bcc, fcc, diamond structure generation
#include "project_a_1.h"
// Project A part 2 functions are reused
// for generating neighbour list for 2 body potentials
#include "project_a_2.h"
using namespace std;

// Argon LJ params: https://www.researchgate.net/figure/Lennard-Jones-LJ-potential-parameters-of-different-materials-considered-in-thepresent_tbl2_319412425
// Argon structure and lattice constant: https://aip.scitation.org/doi/10.1063/1.1726009?cookieSet=1 
// Argon cohesive energy: https://arxiv.org/pdf/2012.05413.pdf 
// https://www.knowledgedoor.com/2/elements_handbook/cohesive_energy.html
//       | epsilon (eV) | sigma (nm) | a (Angstrom) | structure | cohesive energy (eV) |
// Ar-Ar | 0.34         | 0.0104     | 5.311        | FCC       | 0.08                 |
const double pi = 3.141592654;
const double e_0 = 55.26349406 / 1e4; // e^2 eV^-1 Angstrom^-1, vacuum permittivity
const double eV = 1.602176634e-19;

const double a_ArAr = 5.311; // Angstrom, lattice constant
const double epsilon_ArAr = 0.0104; // eV, potential params 
const double sigma_ArAr = 0.34 * 10; // Angstrom, potential params 

const double En_Ar = -0.08; // eV, reference cohesive energy

// MgO Buckingham params: https://www.researchgate.net/figure/Parameters-of-the-pair-interaction-potentials-between-the-ions-in-magnesium-oxide-within_tbl1_226793406
// MgO structural properties: https://materialsproject.org/materials/mp-1265/
// MgO Cohesive energy: https://iopscience.iop.org/article/10.1088/1742-6596/377/1/012067/pdf 
//      | a (Angstrom) | structure          | binding energy (eV) |
// Mg-O | 4.26         | FCC (2 atom basis) | 20.1                |
const double a_MgO = 4.26; // Angstrom, lattice constant

const double A_MgMg = 0; // eV, potential params 
const double B_MgMg = 1; // Angstrom^-1, potential params 
const double C_MgMg = 0; // eV Angstrom^6, potential params 

const double A_MgO = 821.6; // eV, potential params 
const double B_MgO = 1. / 0.3242; // Angstrom^-1, potential params 
const double C_MgO = 0; // eV Angstrom^6, potential params 

const double A_OO = 22764; // eV, potential params 
const double B_OO = 1. / 0.149; // Angstrom^-1, potential params 
const double C_OO = 27.88; // eV Angstrom^6, potential params 

const double q_Mg = 2; // e, Mg 2+ ion charge
const double q_O = -2; // e, O 2- ion charge

const double En_MgO = -20.1; // eV, reference binding energy

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
    // Initialize a empty variable to store total energy of Ar structure
    double Ar_tot_en;
    // Initialize a empty variable to store all Ar atoms real space coordinates, and the real lattice vector
    vector <vector <double> > Ar_cell_xyz, Ar_cell_vec;

    cout << "LJ potential" << endl;
    // Looping for different dimension of simulation cells (NxN)
    for (int n_cells = 5; n_cells < 8; n_cells++) {
        Ar_cell_xyz = get_fcc(n_cells, n_cells, n_cells, a_ArAr); // Part 1 function
        cout << n_cells << " by " << n_cells << " Cell, N atoms: " << Ar_cell_xyz.size() << endl;
        Ar_cell_vec = {{n_cells * a_ArAr, 0, 0}, {0, n_cells * a_ArAr, 0}, {0, 0, n_cells * a_ArAr}}; // Setting the lattice vector to appropiate scale
        vector <neighbour_record> nlist = neighbour_list(Ar_cell_vec, Ar_cell_xyz, 10 * a_ArAr); // Part 2 function, with cutoff = 10 * lattice constant
        Ar_tot_en = 0; // Setting the total Ar structure energy as 0
        for (auto record: nlist) {
            Ar_tot_en += LJ(record.squared_distance, epsilon_ArAr, sigma_ArAr); // For each pair of neighbour, compute their LJ potential and add to total energy
        }
        cout << "Computed Total Cohesive Energy of Ar (eV): " << Ar_tot_en << endl << "Computed Cohesive Energy per Ar atom (eV): " << Ar_tot_en / Ar_cell_xyz.size() << endl;
        cout << "Reference value of Cohesive Energy of Ar per atom (eV): " << En_Ar << endl;
    }
        
    cout << endl << "Coulomb-Buckingham Potential" << endl;
    // Initialize a empty variable to store all MgO atoms real space coordinates, and the real lattice vector
    vector <vector <double> > MgO_cell_vec, MgO_cell_xyz, O_cell_xyz;
    // Looping for different dimension of simulation cells (NxN)
    for (int n_cells = 5; n_cells < 8; n_cells++) {
        // Setting the lattice vector to appropiate scale
        MgO_cell_vec = {{n_cells * a_MgO, 0, 0}, {0, n_cells * a_MgO, 0}, {0, 0, n_cells * a_MgO}};
        // The name is MgO cell firstly initialize Mg atoms, part 1 function
        MgO_cell_xyz = get_fcc(n_cells, n_cells, n_cells, a_MgO);
        // Mg length is useful in determining atom type and charge type by indexing
        int Mg_length = MgO_cell_xyz.size(); 
        // Then initialize O atoms, part 1 function 
        O_cell_xyz = get_fcc(n_cells, n_cells, n_cells, a_MgO, a_MgO/2);
        // Insert O atoms to the end of Mg atoms to make MgO cell
        MgO_cell_xyz.insert(end(MgO_cell_xyz), begin(O_cell_xyz), end(O_cell_xyz));
        // Part 2 function, with cutoff = 10 * lattice constant
        vector <neighbour_record> MgO_nlist = neighbour_list(MgO_cell_vec, MgO_cell_xyz, 10 * a_MgO);
        // Setting the total MgO structure energy as 0
        double MgO_tot_en = 0;
        
        cout << n_cells << " by " << n_cells << " Cell, N atoms: " << MgO_cell_xyz.size() << endl;
        for (auto record: MgO_nlist) {
            // For each pair of neighbour, consider different interaction cases
            // Use different sets of params in CoulombBuckingham potential and add to the total MgO energy

            // Because the way I set up the atom list, labels at [0, Mg_length-1] are Mg atoms, [Mg_length, list_size-1] are O atoms
            // Both ions are Mg
            if ((record.label1 < Mg_length) && (record.label2 < Mg_length)) {
                MgO_tot_en += CoulombBuckingham(record.squared_distance, q_Mg, q_Mg, A_MgMg, B_MgMg, C_MgMg);
            }
            // Both ions are O
            else if ((record.label1 > Mg_length) && (record.label2 > Mg_length)) {
                MgO_tot_en += CoulombBuckingham(record.squared_distance, q_O, q_O, A_OO, B_OO, C_OO);
            }
            // One Mg One O
            else {
                MgO_tot_en += CoulombBuckingham(record.squared_distance, q_Mg, q_O, A_MgO, B_MgO, C_MgO);
            }
        }
        cout << "Computed Total Binding Energy of MgO (eV): " << MgO_tot_en << endl << "Computed Binding Energy per MgO atom (eV): " << MgO_tot_en / MgO_cell_xyz.size() << endl;
        cout << "Reference value of Binding Energy of MgO per atom (eV): " << En_MgO << endl;
    }
}