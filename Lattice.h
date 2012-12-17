#include "tbb/tbb.h"
using namespace tbb;
using namespace std;
#include <stdio.h>
#ifdef OS_WINDOWS
#include "io.h"
#include "sys\stat.h"
#include "glut.h"
#elif defined(OS_LINUX)
#include <sys/stat.h>
#include <GL/glut.h>
#endif
#include <stdlib.h>
#include "fcntl.h"

#include "Sachse_fibroblast.h"
#include "LR_cell.h"
#include <iostream>
#include <time.h>

#include "IsoSurface.hpp"
#include "tbb/blocked_range3d.h"
//#include "advisor-annotate.h"


//This header file defines all the biophysically relevant parameters of the systeml; the system size; the functions carrying out the
//calculations as well as the realizations of the parralel_for classes used in the parallel version of the application
extern tbb::atomic<double> error;
extern int N; //the width and (x and y coordinates)
extern int H; //the height (z coordinate)
extern int Redisplay_time;
const double h = 1e-4; //100mkM - space descritezation step size
const double A_Cap = 1.534e-4; //cm2 - Luo Rudy 1994 - cardiac myocite membrane capacity
const double A_Cap_f = 1.534e-4; //cm2 - Luo Rudy 1994 - cardiac myocite membrane capacity
const int N_fib_myo = 5;// the fibroblast per myocite number ratio
const double Vol_myo_Vol_fib = 59.7; // the ratio between the spaces occupied by fibroblast and myocite
const double mult = 0.5;
const double SM_X = 0.1*mult;//S/m - basal value of myocite tissue electrical conductance in the x direction
const double SM_Z = 0.1*mult; //S/m
const double SM_Y = 0.2*mult; //S/m
const double SF_X = 0.02*mult;// S/m
const double SF_Y = 0.02*mult;// S/m
const double SF_Z = 0.02*mult;// S/m
const double sm_x = SM_X*0.8/(1.+N_fib_myo/Vol_myo_Vol_fib); // effective electrical conductance in the myocite domain in the x direction
const double sm_z = SM_Z*0.8/(1.+N_fib_myo/Vol_myo_Vol_fib);
const double sm_y = SM_Y*0.8/(1.+N_fib_myo/Vol_myo_Vol_fib);
const double sf_x = SF_X*0.8/(1.+Vol_myo_Vol_fib/N_fib_myo);
const double sf_y = SF_Y*0.8/(1.+Vol_myo_Vol_fib/N_fib_myo);
const double sf_z = SF_Z*0.8/(1.+Vol_myo_Vol_fib/N_fib_myo);
const double se_x = 0.6*mult;//
const double se_y = 0.6*mult;//
const double se_z = 0.6*mult;//
const double Vol_myo = 16e-15; //m3 - myocite volume
const double Vol_fib = 0.268e-15; //m3 - fibroblast volume
const double Betta_myo = (0.8/(1.+N_fib_myo/Vol_myo_Vol_fib))*A_Cap/(Vol_myo*1e6);
const double Betta_fib = (0.8/(1.+Vol_myo_Vol_fib/N_fib_myo))*A_Cap_f/(Vol_myo*1e6); //1/m3
const double R_myo_fib = 1e9; //Om myocite-fibroblast interdomain resistance
const double nav = 1; //average number of gap junctions between a fibroblast and a myocite
const double Betta_myo_fib  = 0.8/((Vol_myo+Vol_fib*N_fib_myo))*nav*N_fib_myo;
extern double dt;

void Init(double *Vm, double *Vf, double *mG, double *hG, double *jG, double *dG, double *fG,
          double *XG, double *Cai, double *fs_m, double *fs_f, double *fs_e, Fibroblast *FB, double *f_sum, double *L_Fe);
void OdeSolve_myocyte(double &Vm, double &mG, double &hG, double &jG, double &dG, double &fG, double &XG, double &Cai, double &dVm);
inline int Substeps(double &vd);//devides step length due to value of Voltage (V)
double SolveEquations(double *Vm, double *Vf, double *mG, double *hG, double *jG, double *dG,
                      double *fG, double *XG, double *Cai, double *fs_m, double *fs_f, double *fs_e, double *f_sum,
                      double *L_Vm, double *L_Vf,double *L_Fe1, double *L_Fe2, double *I_me,double *If_e,
                      Fibroblast *FB, double *buff, double *Fe);//Solves the task
void OdeSolve_fib(int i,  double *V, Fibroblast *FB);

//USED IN THE SERIAL VERSION OF THE APPLICATION ONLY
void Get_f_sum(double *Vm, double *Vf, double *fs_m, double *fs_f, double *fs_e, double *f_sum, double *L_Vm, double *L_Vf);
void Get_external_currents(double *Vm, double *Vf, double *Fe, double *fs_m, double *fs_f, double *Im_e,
                           double *If_e, double *L_Vm, double *L_Vf, double *L_Fe1, double *L_Fe2);
void SolvePoisson(double *f_sum, int N, int H, double h, double ssx, double ssy, double ssz, double *buff, double *Fe);
//////////////////////////////////////////////////////


//PARALLEL VERSIONS
class Psolve{
    double *Vm, *mG, *hG, *jG, *dG, *fG, *XG, *Cai, *dVm;
    double *Vf, *dVf;
    Fibroblast *FB;
    double *Im_e, *If_e;
public:
void operator()(const blocked_range<int>& r)const{
    //int n1, n2, n3;//DELETE
        for (int j=r.begin(); j<r.end(); j++)
        {
            OdeSolve_myocyte(Vm[j],mG[j],hG[j],jG[j],dG[j],fG[j],XG[j],Cai[j],dVm[j]);
            dVf[j] = Vf[j];
            OdeSolve_fib(j,Vf,FB);
            dVf[j] = (Vf[j]-dVf[j])/dt+If_e[j];
            Vm[j]+=dt*Im_e[j];
            dVm[j] += Im_e[j];
            Vf[j]+=dt*If_e[j];
        }
    }
    Psolve(double *Vi, double *mGi, double *hGi, double  *jGi, double  *dGi, double *fGi, double *XGi, double *Caii,
        double *Vfi, Fibroblast *FBi, double *Im_ei, double *If_ei, double *dVmi, double *dVfi):Vm(Vi),mG(mGi),hG(hGi),
        jG(jGi),dG(dGi),fG(fGi),XG(XGi),Cai(Caii),Vf(Vfi),FB(FBi),Im_e(Im_ei),If_e(If_ei),dVm(dVmi),dVf(dVfi){}
};

class GetTotalCurrent{
    double *Vm, *Vf, *fs_m, *fs_f, *fs_e, *f_sum, *L_Vm, *L_Vf;
public:
void operator()(const blocked_range<int>& r)const{
    int n1, n2, n3;
    int rn, ln, un, dn, tn, bn;
        for (int j=r.begin(); j<r.end(); j++)
        {
            n3 = j/(N*N);
            n1 = (j-N*N*n3)/N;
            n2 = j-n1*N-N*N*n3;

            if (n1 > 0) ln = n3*N*N + (n1-1)*N + n2;
            else ln = j;

            if (n1 < N-1) rn = n3*N*N + (n1+1)*N + n2;
            else rn = j;

            if (n2 > 0) un = n3*N*N + n1*N + n2 - 1;
            else un = j;

            if (n2 < N-1) dn = n3*N*N + n1*N + n2 + 1;
            else dn = j;

            if (n3 > 0)
                tn = (n3-1)*N*N + n1*N + n2;
            else
                tn = j;

            if (n3 < H-1) bn = (n3+1)*N*N + n1*N + n2;
            else bn = j;
            L_Vm[j] = 1e-3*(sm_x*(Vm[rn]+Vm[ln]-2.*Vm[j])+sm_y*(Vm[un]+Vm[dn]-2.*Vm[j])+sm_z*(Vm[tn]+Vm[bn]-2.*Vm[j]))/(h*h);
            L_Vf[j] = 1e-3*(sf_x*(Vf[rn]+Vf[ln]-2.*Vf[j])+sf_y*(Vf[un]+Vf[dn]-2.*Vf[j])+sf_z*(Vf[tn]+Vf[bn]-2.*Vf[j]))/(h*h);
            f_sum[j] = -L_Vm[j]-L_Vf[j]-fs_m[j]-fs_f[j]-fs_e[j];
        }
    }
    GetTotalCurrent(double *Vmi, double *Vfi, double *fs_mi, double *fs_fi, double *fs_ei, double *f_sumi, double *L_Vmi, double *L_Vfi):Vm(Vmi),
        Vf(Vfi),fs_m(fs_mi),fs_f(fs_fi),fs_e(fs_ei),f_sum(f_sumi),L_Vm(L_Vmi),L_Vf(L_Vfi){
    }
};

class GetExternalCurrents{
      double *Vm,   *Vf,   *Fe,   *fs_m,   *fs_f,   *Im_e,   *If_e,   *L_Vm,   *L_Vf,   *L_Fe1, *L_Fe2;
public:
void operator()(const blocked_range<int>& r)const{
    int n1, n2, n3;
    int rn, ln, un, dn, vn, bn, s;
        for (int j=r.begin(); j<r.end(); j++)
        {
        n3 = j/(N*N);
        n1 = (j-N*N*n3)/N;
        n2 = j-n1*N-N*N*n3;

        n3++;n1++;n2++;

        rn = n3*(N+2)*(N+2) + (n1+1)*(N+2) + n2;
        ln = n3*(N+2)*(N+2) + (n1-1)*(N+2) + n2;
        un = n3*(N+2)*(N+2) + n1*(N+2) + n2 - 1;
        dn = n3*(N+2)*(N+2) + n1*(N+2) + n2 + 1;
        vn = (n3-1)*(N+2)*(N+2) + n1*(N+2) + n2;
        bn = (n3+1)*(N+2)*(N+2) + n1*(N+2) + n2;
        s = n3*(N+2)*(N+2) + n1*(N+2) + n2;

        L_Fe1[j] = 1e-3*(sm_x*(Fe[rn]+Fe[ln]-2.*Fe[s])+sm_y*(Fe[un]+Fe[dn]-2.*Fe[s])+sm_z*(Fe[vn]+Fe[bn]-2.*Fe[s]))/(h*h);
        L_Fe2[j] = 1e-3*(sf_x*(Fe[rn]+Fe[ln]-2.*Fe[s])+sf_y*(Fe[un]+Fe[dn]-2.*Fe[s])+sf_z*(Fe[vn]+Fe[bn]-2.*Fe[s]))/(h*h);
        Im_e[j] = (L_Vm[j]+L_Fe1[j]+fs_m[j]-Betta_myo_fib*(Vm[j]-Vf[j])/(1000.*R_myo_fib))/Betta_myo;
        If_e[j] = (L_Vf[j]+L_Fe2[j]+fs_f[j]+Betta_myo_fib*(Vm[j]-Vf[j])/(1000.*R_myo_fib))/Betta_fib;
        }
    }
    GetExternalCurrents(double *Vmi, double *Vfi, double *Fei, double *fs_mi, double *fs_fi, double *Im_ei, double *If_ei,
        double *L_Vmi, double *L_Vfi, double *L_Fe1i, double *L_Fe2i):
     Vm(Vmi),Vf(Vfi),Fe(Fei),fs_m(fs_mi),fs_f(fs_fi),Im_e(Im_ei),If_e(If_ei),L_Vm(L_Vmi),L_Vf(L_Vfi),L_Fe1(L_Fe1i),L_Fe2(L_Fe2i){
    }
};

class PoissonSolver{
    double *f_sum,   ssx,   ssy,   ssz,   *buff,   *Fe;
     //tbb::atomic<double*> error;

public:
void operator()(const blocked_range<int>& r)const{
    tbb::atomic<double>  er;
    int n1,n2,n3;
    int rn, ln, un, dn, vn, bn;
    tbb::atomic<double> error_local;
    error_local = 0;
        for (int j=r.begin(); j<r.end(); j++)
        {
        n3 = j/(N*N);
        n1 = (j-N*N*n3)/N;
        n2 = j-n1*N-N*N*n3;

        n3++;n1++;n2++;

        rn = n3*(N+2)*(N+2) + (n1+1)*(N+2) + n2;
        ln = n3*(N+2)*(N+2) + (n1-1)*(N+2) + n2;
        un = n3*(N+2)*(N+2) + n1*(N+2) + n2 - 1;
        dn = n3*(N+2)*(N+2) + n1*(N+2) + n2 + 1;
        vn = (n3-1)*(N+2)*(N+2) + n1*(N+2) + n2;
        bn = (n3+1)*(N+2)*(N+2) + n1*(N+2) + n2;
        buff[j] = (ssz*(Fe[vn]+Fe[bn])+ssx*(Fe[rn]+Fe[ln])+ssy*(Fe[un]+Fe[dn])-f_sum[j]*h*h)/(2*(ssx+ssy+ssz));

        er = fabs(buff[j]-Fe[n3*(N+2)*(N+2)+n1*(N+2)+n2]);
        if (er > error_local) error_local = er;
        }

        if (error_local > error) error = error_local;
    }
    //PoissonSolver(double *f_sumi, double ssxi, double ssyi, double sszi, double *buffi, double *Fei, tbb::atomic<double*> errori):
     //f_sum(f_sumi),ssx(ssxi),ssy(ssyi),ssz(sszi),buff(buffi),Fe(Fei),error(errori){}
PoissonSolver(double *f_sumi, double ssxi, double ssyi, double sszi, double *buffi, double *Fei):
     f_sum(f_sumi),ssx(ssxi),ssy(ssyi),ssz(sszi),buff(buffi),Fe(Fei){}
};

class SolvePoissonUpdate{
      double *buff,   *Fe;;
public:
void operator()(const blocked_range<int>& r)const{

    int n1,n2,n3;
        for (int j=r.begin(); j<r.end(); j++)
        {
        n3 = j/(N*N);
        n1 = (j-N*N*n3)/N;
        n2 = j-n1*N-N*N*n3;
        n3++;n1++;n2++;
        Fe[n3*(N+2)*(N+2)+n1*(N+2)+n2] = buff[j];
        }
    }
    SolvePoissonUpdate(double *buffi, double *Fei):
     buff(buffi),Fe(Fei){
    }
};

class getVoxels
{
    double *Vm;
    double ***Voxels;
public:
    void operator()(const blocked_range3d<int>& r)const{
        for (int i = r.pages().begin(); i<r.pages().end(); i++)
            for (int j = r.rows().begin(); j<r.rows().end(); j++)
                for (int k = r.cols().begin(); k<r.cols().end(); k++)
                    Voxels[i][j][k] = Vm[k*N*N+j*N+i];
    }
    getVoxels(double *Vmi, double ***Voxelsi):Vm(Vmi),Voxels(Voxelsi){}
};

class findIsoSurface
{
    double ***voxels;
    double *Der;
    double isovalue;
    vector<vertex> vertexList;
    int sign;
public:
    void operator()(const blocked_range3d<int>& r)const{
            for (int x = r.pages().begin(); x<r.pages().end(); x++)
                for (int y = r.rows().begin(); y<r.rows().end(); y++)
                    for (int z = r.cols().begin(); z<r.cols().end(); z++)
                    {

                    }
    }
    findIsoSurface(double ***voxelsi, double *Deri, double isovaluei):
      voxels(voxelsi), Der(Deri), isovalue(isovaluei){}
};
extern vector<vertex> vertices_m;
extern vector<vertex> vertices_f;
extern bool DrawReady;
extern double speed;//time itarator
extern double ***voxels_m;
extern double ***voxels_f;
extern double *dVm;
extern double *dVf;
