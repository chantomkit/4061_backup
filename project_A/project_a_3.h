#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;

// Argon LJ params: https://www.researchgate.net/figure/Lennard-Jones-LJ-potential-parameters-of-different-materials-considered-in-thepresent_tbl2_319412425
// Argon structure and lattice constant: https://aip.scitation.org/doi/10.1063/1.1726009?cookieSet=1 
// Argon cohesive energy: https://www.knowledgedoor.com/2/elements_handbook/cohesive_energy.html
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
// MgO binding energy: https://iopscience.iop.org/article/10.1088/1742-6596/377/1/012067/pdf 
//      | a (Angstrom) | structure          | binding energy (eV) |
// Mg-O | 4.26         | FCC (2 atom basis) | 20.1                |
const double a_MgO = 4.26; // Angstrom, lattice constant
// const double a_MgO = 1; // Angstrom, lattice constant

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