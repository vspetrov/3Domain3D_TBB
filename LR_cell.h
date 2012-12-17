#ifndef LR_cell_H
#define LR_cell_H
#include "math.h"


//right-hand parts of the Luo-Rudy System
extern double VFunction(double V, double mG, double hG, double jG, double dG, double fG, double XG, double Cai);
extern double mFunction(double V, double mG, double &delta_t);
extern double hFunction(double V, double hG, double &delta_t);
extern double jFunction(double V, double jG, double &delta_t);
extern double dFunction(double V, double dG, double &delta_t);
extern double fFunction(double V, double fG, double &delta_t);
extern double XFunction(double V, double XG, double &delta_t);
extern double CaiFunction(double Caii,double dG, double fG, double V);

#endif