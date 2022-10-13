#include <iostream>
#include <math.h>
#include <vector>
using namespace std;

double tr(vector <vector <double>> a) {
    double sum = 0;
    for (int i = 0; i < a.size(); i++) {
        sum += a[i][i];
    }
    return sum;
}

vector <vector <double>> matProd(vector <vector <double>> a, vector <vector <double>> b) {
    vector <vector <double>> c(a.size(), vector <double>(a.size(), 0));
    for (int i = 0; i < a.size(); i++) {
        for (int j = 0; j < a.size(); j++) {
            for (int k = 0; k < a.size(); k++)
                c[i][j] += a[i][k] * b[k][j];
        }
    }
    return c;
}

vector <vector <double>> matInverse(vector <vector <double>> a) {
    int n = a.size();
    vector <double> c(n, 0);
    vector <vector <vector <double>>> s(n, vector <vector <double>>(n, vector <double>(n, 0)));
    vector <vector <double>> d(n, vector <double>(n, 0));

    for (int i = 0; i < n; i++)
        s[0][i][i] = 1;

    for (int k = 1; k < n; k++ ) {
        s[k] = matProd(a, s[k - 1]);
        c[n - k] = -tr(s[k]) / k;
        for (int i = 0; i < n; i++)
            s[k][i][i] += c[n-k];
    }

    c[0] = -tr(matProd(a, s[n-1])) / n;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            d[i][j] = -s[n - 1][i][j] / c[0];
    return d;
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

int main() {
    vector <vector <double>> a {{1,2,3,4},{4,3,2,1},{1,4,2,3},{2,1,3,4}};
    int n = a.size();
    vector <vector <double>> d = matInverse(a);
    cout << "Original Matrix" << endl;
    print_2dvector(a);
    cout << "Matrix Inverse" << endl;
    print_2dvector(d);
}