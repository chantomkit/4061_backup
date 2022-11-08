#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define n 100
#define j 2


#define g 9.80
#define speed 67
#define mass 250
#define area 0.93
#define density 1.2

void jump_parallelworld(double);
void jump_originalworld(double);
double deg2rad(double);

int main() {
    jump_parallelworld(deg2rad(42.5));
    // jump_parallelworld(deg2rad(45.5));
    // jump_parallelworld(deg2rad(47.5));

    // jump_originalworld(deg2rad(42.5));
    // jump_originalworld(deg2rad(45.5));
    // jump_originalworld(deg2rad(47.5));
}

double deg2rad(double angle) {
    return angle * M_PI / 180;
}

void jump_parallelworld(double angle) {
    static double x[n+1], y[n+1], vx[n+1], vy[n+1], ax[n+1], ay[n+1];
    double k = area * density / (2 * mass);
    double dt = 2 * speed * sin(angle) / (g * n);
    double d = dt * dt / 2;

    x[0] = y[0] = 0;
    vx[0] = speed * cos(angle);
    vy[0] = speed * sin(angle);
    double v = sqrt(vx[0] * vx[0] + vy[0] * vy[0]);
    ax[0] = -k * pow(v, 2/5) * vx[0];
    ay[0] = -g - k * pow(v, 2/5) * vy[0];

    double p = vx[0] * ax[0] + vy[0] * ay[0];
    x[1] = x[0] + dt * vx[0] + d * ax[0];
    y[1] = y[0] + dt * vy[0] + d * ay[0];
    vx[1] = vx[0] + dt * ax[0] - d * k * (pow(v, 2/5) * ax[0] + 2 * p * vx[0] / (5 * pow(v, 8/5)));
    vy[1] = vy[0] + dt * ay[0] - d * k * (pow(v, 2/5) * ay[0] + 2 * p * vy[0] / (5 * pow(v, 8/5)));

    v = sqrt(vx[1] * vx[1] + vy[1] * vy[1]);
    ax[1] = -k * v * vx[1];
    ay[1] = -g - k * v * vy[1];

    double d2 = 2 * dt;
    double d3 = dt / 3;

    for (int i = 0; i <= (n - 2); i++) {
        x[i+2] = x[i] + d2 * vx[i+1];
        y[i+2] = y[i] + d2 * vy[i+1];
        vx[i+2] = vx[i] + d2 * ax[i+1];
        vy[i+2] = vy[i] + d2 * ay[i+1];
        v = sqrt(vx[i+2] * vx[i+2] + vy[i+2] * vy[i+2]);
        ax[i+2] = -k * pow(v, 2/5) * vx[i+2];
        ay[i+2] = -g - k * pow(v, 2/5) * vy[i+2];

        x[i+2] = x[i] + d3 * (vx[i+2] + 4 * vx[i+1] + vx[i]);
        y[i+2] = y[i] + d3 * (vy[i+2] + 4 * vy[i+1] + vy[i]);
        vx[i+2] = vx[i] + d3 * (ax[i+2] + 4 * ax[i+1] + ax[i]);
        vy[i+2] = vy[i] + d3 * (ay[i+2] + 4 * ay[i+1] + ay[i]);
    }

    for (int i = 0; i <= n; i++) {
        printf("%f %f\n", x[i], y[i]);
    }
    return;
}

void jump_originalworld(double angle) {
    static double x[n+1], y[n+1], vx[n+1], vy[n+1], ax[n+1], ay[n+1];
    double k = area * density / (2 * mass);
    double dt = 2 * speed * sin(angle) / (g * n);
    double d = dt * dt / 2;

    x[0] = y[0] = 0;
    vx[0] = speed * cos(angle);
    vy[0] = speed * sin(angle);
    double v = sqrt(vx[0] * vx[0] + vy[0] * vy[0]);
    ax[0] = -k * v * vx[0];
    ay[0] = -g - k * v * vy[0];

    double p = vx[0] * ax[0] + vy[0] * ay[0];
    x[1] = x[0] + dt * vx[0] + d * ax[0];
    y[1] = y[0] + dt * vy[0] + d * ay[0];
    vx[1] = vx[0] + dt * ax[0] - d * k * (v * ax[0] + p * vx[0] / v);
    vy[1] = vy[0] + dt * ay[0] - d * k * (v * ay[0] + p * vy[0] / v);

    v = sqrt(vx[1] * vx[1] + vy[1] * vy[1]);
    ax[1] = -k * v * vx[1];
    ay[1] = -g - k * v * vy[1];

    double d2 = 2 * dt;
    double d3 = dt / 3;

    for (int i = 0; i <= (n - 2); i++) {
        x[i+2] = x[i] + d2 * vx[i+1];
        y[i+2] = y[i] + d2 * vy[i+1];
        vx[i+2] = vx[i] + d2 * ax[i+1];
        vy[i+2] = vy[i] + d2 * ay[i+1];
        v = sqrt(vx[i+2] * vx[i+2] + vy[i+2] * vy[i+2]);
        ax[i+2] = -k * v * vx[i+2];
        ay[i+2] = -g - k * v * vy[i+2];

        x[i+2] = x[i] + d3 * (vx[i+2] + 4 * vx[i+1] + vx[i]);
        y[i+2] = y[i] + d3 * (vy[i+2] + 4 * vy[i+1] + vy[i]);
        vx[i+2] = vx[i] + d3 * (ax[i+2] + 4 * ax[i+1] + ax[i]);
        vy[i+2] = vy[i] + d3 * (ay[i+2] + 4 * ay[i+1] + ay[i]);
    }

    for (int i = 0; i <= n; i++) {
        printf("%f %f\n", x[i], y[i]);
    }
    return;
}