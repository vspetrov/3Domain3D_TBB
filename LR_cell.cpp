#include "LR_cell.h"
//in this file there are the realizations of the right hand side functions that difine the dynamics of cardiac MYOCITE
// according to the model of Sachse et. al. These functions are used in the routine 
//"void OdeSolve_myocyte(double &Vm, double &mG, double &hG, double &jG, double &dG, double &fG, double &XG, double &Cai)"
//integrating the myocite dynamics.

//cardiac myocite biophysical cell constants
const double G_Na=22.; //maximum conductance of the sodium ionic channel
const double E_Na=54.7939; //Nernst potential for the sodium ionic channel
const double G_K=1; //maximum conductance of the stationary potassium ionic channel
const double E_Kl=-85.2268;//Nernst potential for the non-stationary potassium ionic channel
const double E_k=-77.58;//Nernst potential for the stationary potassium ionic channel
const double G_Kl=0.6047;//maximum conductance of the non-stationary potassium ionic channel
const double G_si=0.02;//maximum conductance of the slow inward calcium ionic channel
const double G_b = 0.03921;//maximum conductance of the leakage

double VFunction(double V, double mG, double hG, double jG, double dG, double fG, double XG, double Cai) //find V for the Cell numbered CN_c, CN_r
{
	double i_Na;
	double i_si;
	double E_si;
	double i_K;
	double Xi;
	double i_Kl;
	double Kl_inf;
	double alpha_Kl, betta_Kl;
	double i_Kp;
	double Kp;
	double i_b;
	//find fast sodium current I_Na
	i_Na=G_Na*mG*mG*mG*hG*jG*(V-E_Na);
	
	//find slow-inward current I_si
	E_si=7.7-13.0287*log(Cai);
	i_si=G_si*dG*fG*(V-E_si);

	//find time-dependent potassium current I_K
	if (V>-100.)
		Xi=2.837*(exp(0.04*(V+77.))-1.)/((V+77.)*exp(0.04*(V+35.)));
	else
		Xi=1.;
	i_K=G_K*XG*Xi*(V-E_k);
	
	//find time-independent potassium current I_K1
	alpha_Kl=1.02/(1+exp(0.2385*(V-E_Kl-59.215)));
	betta_Kl=(0.49124*exp(0.08032*(V-E_Kl+5.476))+
			exp(0.06175*(V-E_Kl-594.31)))/(1.+exp(-0.5143*(V-E_Kl+4.753)));
	Kl_inf=alpha_Kl/(alpha_Kl+betta_Kl);
	i_Kl=G_Kl*Kl_inf*(V-E_Kl);
		
	//find platou potassium current I_Kp
	Kp=1./(1.+exp((7.488-V)/5.98));
	i_Kp=0.0183*Kp*(V-E_Kl);
		
	//find background current I_b
	
	i_b=G_b*(V+59.87);

	return -(i_Na+i_si+i_K+i_Kl+i_Kp+i_b);
}

double mFunction(double V, double mG, double &delta_t)
{
	double p,q;
	double alpha, betta;
	alpha=0.32*(V+47.13)/(1.-exp(-0.1*(V+47.13)));
	betta=0.08*exp(-V/11.);
	p=alpha/(alpha+betta);
	q=1./(alpha+betta);
	return p+(mG-p)*exp(-delta_t/q);
}

double hFunction(double V, double hG, double &delta_t)
{
	double p,q;
	double alpha, betta;

	if (V < -40.) 
		{
			alpha=0.135*exp((80.+V)/-6.8);
			betta=3.56*exp(0.079*V)+3.1e5*exp(0.35*V);
		}
	else
		{
			alpha=0.;
			betta=1./(0.13*(1.+exp((V+10.66)/-11.1)));
		}
	p=alpha/(alpha+betta);
	q=1./(alpha+betta);
	return p+(hG-p)*exp(-delta_t/q);
}

double jFunction(double V, double jG, double &delta_t)
{
	double p,q;
	double alpha, betta;
	if (V < -40.) 
		{
			alpha=(-1.2714e5*exp(0.2444*V)-3.474e-5*exp(-0.04391*V))*
			(V+37.78)/(1+exp(0.311*(V+79.23)));
			betta=0.1212*exp(-0.01052*V)/(1.+exp(-0.1378*(V+40.14)));
		}
	else
		{
			alpha=0.;
			betta=0.3*exp(-2.535e-7*V)/(1.+exp(-0.1*(V+32.)));
		}
	p=alpha/(alpha+betta);
	q=1./(alpha+betta);
	return p+(jG-p)*exp(-delta_t/q);
}

double dFunction(double V, double dG, double &delta_t)
{
	double p,q;
	double alpha, betta;
	alpha=0.095*exp(-0.01*(V-5.0))/(1.+exp(-0.072*(V-5.0)));
	betta=0.07*exp(-0.017*(V+44.0))/(1.+exp(0.05*(V+44.0)));
	p=alpha/(alpha+betta);
	q=1./(alpha+betta);
	return p+(dG-p)*exp(-delta_t/q);
}

double fFunction(double V, double fG, double &delta_t)
{
	double p,q;
	double alpha, betta;
	alpha=0.012*exp(-0.008*(V+28.))/(1.+exp(0.15*(V+28.)));
	betta=0.0065*exp(-0.02*(V+30.))/(1.+exp(-0.2*(V+30.)));
	p=alpha/(alpha+betta);
	q=1./(alpha+betta);
	return p+(fG-p)*exp(-delta_t/q);
}

double XFunction(double V, double XG, double &delta_t)
{
	double p,q;
	double alpha, betta;
	alpha=0.0005*exp(0.083*(V+50.))/(1.+exp(0.057*(V+50.)));
	betta=0.0013*exp(-0.06*(V+20.))/(1.+exp(-0.04*(V+20.)));
	p=alpha/(alpha+betta);
	q=1./(alpha+betta);
	return p+(XG-p)*exp(-delta_t/q);
}

double CaiFunction(double Cai, double dG, double fG, double V)
{
	double i_si;
	double E_si;
	E_si=7.7-13.0287*log(Cai);
	i_si=G_si*dG*fG*(V-E_si);
	return -1e-4*i_si+0.07*(1e-4-Cai);
}