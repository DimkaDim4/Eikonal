#include <iostream>
#include "Solver.h"

int I = 500;
double h = 1. / (I - 1.);

int main(int argc, char* argv[])
{
    double* vp = new double[I * I * I];
    double* vs = new double[I * I * I];
    double* rho = new double[I * I * I];

    for (int i = 0; i < I; i++)
    {
        for (int j = 0; j < I; j++)
        {
            for (int k = 0; k < I / 3; k++)
            {
                vp[(i * I + j) * I + k] = 3.2;
                vs[(i * I + j) * I + k] = 1.82;
                rho[(i * I + j) * I + k] = 2.7;
            }
            for (int k = I / 3; k < 2 * I / 3; k++)
            {
                vp[(i * I + j) * I + k] = 5.9;
                vs[(i * I + j) * I + k] = 3.42;
                rho[(i * I + j) * I + k] = 2.85;
            }
            for (int k = 2 * I / 3; k < I; k++)
            {
                vp[(i * I + j) * I + k] = 6.95;
                vs[(i * I + j) * I + k] = 4.03;
                rho[(i * I + j) * I + k] = 2.81;
            }
        }
    }

    /*for (int i = 0; i < I; i++)
    {
        for (int j = 0; j < I; j++)
        {
            for (int k = 0; k < I; k++)
            {
                if (k < I / 2)
                {
                    vp[(i * I + j) * I + k] = 1;
                    vs[(i * I + j) * I + k] = 1;
                    rho[(i * I + j) * I + k] = 1;
                }
                else
                {
                    if (i < I / 2)
                    {
                        if (j < I / 2)
                        {
                            vp[(i * I + j) * I + k] = 2;
                            vs[(i * I + j) * I + k] = 2;
                            rho[(i * I + j) * I + k] = 2;
                        }
                        else
                        {
                            vp[(i * I + j) * I + k] = 3;
                            vs[(i * I + j) * I + k] = 3;
                            rho[(i * I + j) * I + k] = 3;
                        }
                    }
                    else
                    {
                        if (j < I / 2)
                        {
                            vp[(i * I + j) * I + k] = 4;
                            vs[(i * I + j) * I + k] = 4;
                            rho[(i * I + j) * I + k] = 4;
                        }
                        else
                        {
                            vp[(i * I + j) * I + k] = 5;
                            vs[(i * I + j) * I + k] = 5;
                            rho[(i * I + j) * I + k] = 5;
                        }
                    }
                }


            }
        }
    }*/

    //std::vector<std::vector<bool>> g = gen(4);

    EnviromentModel model;
    model.InitSize(I, I, I);
    model.InitEnviromentFromArray(vp, vs, rho);

    const std::list<ContactBoundary>* cb = model.Bounds();

    Eikonal eik;
    eik.SetModel(&model);
    eik.AddSourse(0.5, 0.5, 0.2);
    eik.Calculate(0.6);
    

    //Wave3d wave(vp, vs, rho, I, 1.2, f);

    delete[] vp;
    delete[] vs;
    delete[] rho;

    //wave.solve();
}
