#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <limits>
// Project A part 1 functions are reused
// simple cubic, bcc, fcc, diamond structure generation
#include "project_a_1.h"
// Project A part 2 functions are reused
// for generating neighbour list for 2 body potentials
#include "project_a_2.h"
// Constants and potentials
#include "project_a_3.h"
using namespace std;

int main() {
    // Initialize a empty variable to store total energy of Ar structure
    double Ar_tot_en;
    // Initialize a empty variable to store all Ar atoms real space coordinates, and the real lattice vector
    vector <vector <double> > Ar_cell_xyz, Ar_cell_vec;

    cout << "LJ potential" << endl;
    // Looping for different dimension of simulation cells (NxN)
    for (int dim = 1; dim < 5; dim++) {
        int n_cells = dim*2;
        Ar_cell_xyz = get_fcc(n_cells, n_cells, n_cells, a_ArAr); // Part 1 function
        cout << n_cells << " by " << n_cells << " by " << n_cells << " Cell, N atoms: " << Ar_cell_xyz.size() << endl;
        Ar_cell_vec = {{n_cells * a_ArAr, 0, 0}, {0, n_cells * a_ArAr, 0}, {0, 0, n_cells * a_ArAr}}; // Setting the lattice vector to appropiate scale
        vector <neighbour_record> nlist = neighbour_list(Ar_cell_vec, Ar_cell_xyz, numeric_limits<int>::max()); // Part 2 function, with cutoff = 10 * lattice constant
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
    for (int dim = 1; dim < 5; dim++) {
        int n_cells = dim*2;
        // Setting the lattice vector to appropiate scale
        MgO_cell_vec = {{n_cells * a_MgO, 0, 0}, {0, n_cells * a_MgO, 0}, {0, 0, n_cells * a_MgO}};
        // The name is MgO cell firstly initialize Mg atoms, part 1 function
        MgO_cell_xyz = get_fcc(n_cells, n_cells, n_cells, a_MgO);
        // Mg length is useful in determining atom type and charge type by indexing
        int Mg_length = MgO_cell_xyz.size(); 
        // Then initialize O atoms, part 1 function 
        O_cell_xyz = get_fcc(n_cells, n_cells, n_cells, a_MgO, 0.5);
        // Insert O atoms to the end of Mg atoms to make MgO cell
        MgO_cell_xyz.insert(end(MgO_cell_xyz), begin(O_cell_xyz), end(O_cell_xyz));
        // Part 2 function, with cutoff = 10 * lattice constant
        vector <neighbour_record> MgO_nlist = neighbour_list(MgO_cell_vec, MgO_cell_xyz, numeric_limits<int>::max());
        // Setting the total MgO structure energy as 0
        double MgO_tot_en = 0;
        
        cout << n_cells << " by " << n_cells <<" by " << n_cells << " Cell, N atoms: " << MgO_cell_xyz.size() << endl;
        for (auto record: MgO_nlist) {
            // For each pair of neighbour, consider different interaction cases
            // Use different sets of params in CoulombBuckingham potential and add to the total MgO energy

            // Because the way I set up the atom list, labels at [0, Mg_length-1] are Mg atoms, [Mg_length, list_size-1] are O atoms
            // Both ions are Mg
            if ((record.label1 < Mg_length) && (record.label2 < Mg_length)) {
                MgO_tot_en += CoulombBuckingham(record.squared_distance, q_Mg, q_Mg, A_MgMg, B_MgMg, C_MgMg);
            }
            // Both ions are O
            else if ((record.label1 >= Mg_length) && (record.label2 >= Mg_length)) {
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