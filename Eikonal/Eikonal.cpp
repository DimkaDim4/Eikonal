#include <iostream>
#include "Solver.h"

int I = 150;
double h = 1. / (I - 1.);

int main(int argc, char* argv[])
{
    double* vp = new double[I * I * I];
    double* vs = new double[I * I * I];
    double* rho = new double[I * I * I];

    Index sphere_1 = Index(0, 0, 0);
    Index sphere_2 = Index(75, 75, 75);

    for (int i = 0; i < I; i++)
    {
        for (int j = 0; j < I; j++)
        {
            /*for (int k = 0; k < I; k++)
            {
                if ((i - sphere_1.i) * (i - sphere_1.i) + (j - sphere_1.j) * (j - sphere_1.j) + (k - sphere_1.k) * (k - sphere_1.k) <= 25) {
                    vp[(i * I + j) * I + k] = 3.2;
                    vs[(i * I + j) * I + k] = 1.82;
                    rho[(i * I + j) * I + k] = 1.7;
                } else if ((i - sphere_2.i) * (i - sphere_2.i) + (j - sphere_2.j) * (j - sphere_2.j) + (k - sphere_2.k) * (k - sphere_2.k) <= 25) {
                    vp[(i * I + j) * I + k] = 0;
                    vs[(i * I + j) * I + k] = 0;
                    rho[(i * I + j) * I + k] = 0;
                } else {
                    vp[(i * I + j) * I + k] = 6.9;
                    vs[(i * I + j) * I + k] = 3.82;
                    rho[(i * I + j) * I + k] = 2.85;
                }
            }*/

            for (int k = 0; k < I / 3; k++)
            {
                vp[(i * I + j) * I + k] = 3.2;
                vs[(i * I + j) * I + k] = 1.82;
                rho[(i * I + j) * I + k] = 1.7;
            }
            for (int k = I / 3; k < 2 * I / 3; k++)
            {
                vp[(i * I + j) * I + k] = 6.9;
                vs[(i * I + j) * I + k] = 3.82;
                rho[(i * I + j) * I + k] = 2.85;
            }
            for (int k = 2 * I / 3; k < I; k++)
            {
                vp[(i * I + j) * I + k] = 2. * 6.9;
                vs[(i * I + j) * I + k] = 2. * 3.82;
                rho[(i * I + j) * I + k] = 2. * 2.85;
            }
        }
    }

    EnviromentModel model;
    model.InitSize(I, I, I);
    model.InitEnviromentFromArray(vp, vs, rho);

    const std::list<ContactBoundary>* cb = model.Bounds();

    Eikonal eik;
    eik.SetModel(&model);
    eik.SetSourse(0, 0, 0);
    eik.Calculate(3.);
    

    //Wave3d wave(vp, vs, rho, I, 1.2, f);

    delete[] vp;
    delete[] vs;
    delete[] rho;

    //wave.solve();
}
