#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>

using namespace std;

const int n = 99;
const int m = 2;

void relax(vector <double> &u, vector <double> &d, vector <double> &s, double p, double del, int nmax) {
    int n = u.size() - 1;
    double q = 1 - p, fi = 0;
    double du = 2 * del;
    int k = 0;
    while ((du > del) && (k < nmax)) {
        du = 0;
        for (int i = 1; i < n; i++) {
            fi = u[i];
            u[i] = p * u[i] + q * ((d[i+1] + d[i]) * u[i+1] + (d[i] + d[i-1]) * u[i-1] + 2 * s[i]) / (4 * d[i]);
            fi = u[i] - fi;
            du += fi * fi;
        }
        du = sqrt(du / n);
        k++;
    }
    if (k == nmax) {
        cout << "Convergence not found after " << nmax << " iterations" << endl;
    }
}

int main() {
    vector <double> u(n+1), d(n+1), s(n+1);
    double l = 3, l2 = l/2, h = l/n, h2 = h*h;
    double x0 = 0.25, x2 = x0 * x0, e0 = 1 / M_E;
    double x = 0, rho = 3, g = 9.8, f0 = 200;
    double y = 1e9 * pow(0.03, 3) * 0.2 / 3;
    double u0 = 0.032, p = 1.5, del = 1e-3;
    int nmax = 100;

    for (int i = 0; i <= n; i++) {
        s[i] = rho * g;
        x = h * i - l2;
        if (fabs(x) < x0) {
            s[i] += f0 * (exp(-x * x / x2) - e0);
        }
        s[i] *= h2 / y;
    }

    for (int i = 1; i < n; i++) {
        x = M_PI * h * i / l;
        u[i] = u0 * sin(x);
        d[i] = 1;
    }
    d[0] = d[n] = 1;
    
    relax(u, d, s , p, del, nmax);

    x = 0;
    double mh = m * h;
    for (int i = 0; i < n; i+=m) {
        cout << x << " " << 100*u[i] << endl;
        x += mh;
    }
}