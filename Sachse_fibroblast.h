#ifndef _FIBROBLAST_H
#define _FIBROBLAST_H
#include "math.h"
//there are the definitions of the cardiac fibroblast state variables; constants and righthand side functions describing dynamics

//cell constants
const double Cf = 4.5;
const double R = 8.31;
const double F = 9.65e1;
const double T = 295.;
const double K_0 = 5.;
const double K_i = 140.;
const double E_k = R*T/F*log(K_0/K_i);
const double G_kir = 1.02;
const double a_kir = 0.94;
const double b_kir = 1.27;
const double kv0 = 30e-3;
const double zv = 1.28;
const double k_v0 = 2e-3;
const double z_v = -1.53;
const double ko = 77e-3;
const double k_o = 18e-3;
const double P_shkr = 5.4e-3;
const double G_b = 6.9e-3;
const double E_b = 0.;
const double FbyRT = F/(R*T);

//the structure describing the fibroblast state variables
struct Fibroblast{
	double  O_shkr;
	double  C0;
	double  C1;
	double  C2;
	double  C3;
	double  C4;
};

//right hand side dynamical functions
double Vf_function(double V, double O_shkr);
double C0_function(double C0, double C1, double V);
double C1_function(double C0, double C1, double C2, double V);
double C2_function(double C1, double C2, double C3, double V);
double C3_function(double C2, double C3, double C4, double V);
double C4_function(double C3, double C4, double O_shkr, double V);
double O_function(double C4, double O_shkr);
#endif