#include "Lattice.h"
vector<vertex> vertices_m;
vector<vertex> vertices_f;
tbb::atomic<double> error;
int N;
int H;
int Redisplay_time;
double speed = 0;
double dt=0.02;
double ***voxels_m;
double ***voxels_f;
double *dVm;
double *dVf;
void Excite(double *Vm)
{
	int xx = floor(1+(double)rand()/(double)RAND_MAX*(N-2));
	int yy = floor(1+(double)rand()/(double)RAND_MAX*(N-2));
	int zz = floor(1+(double)rand()/(double)RAND_MAX*(H-2));
	int spot_radius = 2;
	int xs,ys,zs,xe,ye,ze;
	xs = xx - spot_radius;
	ys = yy - spot_radius;
	zs = zz - spot_radius;
	xe = xx + spot_radius;
	ye = yy + spot_radius;
	ze = zz + spot_radius;

	if (xs < 0) xs = 0;
	if (xe > N-1) xe = N-1;

	if (ys < 0) ys = 0;
	if (ye > N-1) ye = N-1;

	if (zs < 0) zs = 0;
	if (ze > H-1) ze = H-1;


	for (int x = xs; x<=xe; x++)
		for (int y = ys; y<=ye; y++)
			for (int z = zs; z<=ze; z++)
				Vm[z*N*N+y*N+x] = -20;

}
void Init(double *Vm, double *Vf, double *mG, double *hG, double *jG, double *dG, double *fG, 
		  double *XG, double *Cai, double *fs_m, double *fs_f, double *fs_e, Fibroblast *FB, double *f_sum, double *Fe)
{
	//initializing function: sets the initial values of the state variables of the model
	int i;
	srand(int(time(NULL)));  
	for (i=0; i<N*N*H; i++)
	{
			Vm[i] = -80.;
			Vf[i] = -60.;
			fs_m[i] = 0;
			fs_f[i]=0.;
			fs_e[i]=0.;
			mG[i]  = 0.00231609 ;
			hG[i]  = 0.973114   ;
			jG[i]  = 0.84991    ;
			dG[i]  = 0.00434296 ;
			fG[i]  = 0.880756   ;
			XG[i]  = 0.018826   ;
			Cai[i] = 0.000445703; 
			FB[i].C0=0.176258;
			FB[i].C1=0.367597;
			FB[i].C2=0.287479;
			FB[i].C3=0.0998994;
			FB[i].C4=0.0129962;
			FB[i].O_shkr = 0.0555375;
	}


	/*for (i =0; i<3; i++)
		for (int j=N-1; j>N-4; j--)
			for (int k=H-1; k>H-4; k--)
				Vm[k*N*N+j*N+i] = -31;*/

	Excite(Vm);

	for (i=0; i<(N+2)*(N+2)*(H+2); i++)
		Fe[i] = 0;

	/* read the initial conditions from binary
	int fd = open("rst.bin",O_RDWR|O_CREAT | O_BINARY,S_IREAD|S_IWRITE);
	read(fd,Vm,N*N*H*sizeof(double));
	read(fd,mG,N*N*H*sizeof(double));
	read(fd,hG,N*N*H*sizeof(double));
	read(fd,jG,N*N*H*sizeof(double));
	read(fd,dG,N*N*H*sizeof(double));
	read(fd,fG,N*N*H*sizeof(double));
	read(fd,XG,N*N*H*sizeof(double));
	read(fd,Cai,N*N*H*sizeof(double));
	read(fd,Vf,N*N*H*sizeof(double));
	close(fd);*/


}

double SolveEquations(double *Vm, double *Vf,  double *mG, double *hG, double *jG,
					  double *dG, double *fG, double *XG,double *Cai, double *fs_m, double *fs_f, double *fs_e,
					  double *f_sum, double *L_Vm, double *L_Vf,double *L_Fe1, double *L_Fe2,double *Im_e,double *If_e,
					  Fibroblast *FB, double *buff, double *Fe)
{
	//the function solves the system for a given time MaxTime: calls all neseccary functions
	int i; //counting variables
	int Time;//time itarator

	time_t t1,t2;


	//tbb::atomic<double*> error; //error estimatoin parameter for poisson solver



	
	//cleaning out the auxilary arrays
	for (i=0; i<N*N*H; i++)
	{
		Im_e[i] = If_e[i] = f_sum[i] = L_Fe1[i] = L_Fe2[i] = L_Vm[i] = L_Vf[i] = 0;
	}
	
		
	

	bool DrawEvent = false;
	dVm = new double[N*N*H];
	dVf = new double[N*N*H];

	voxels_m = new double**[N];
	voxels_f = new double**[N];
	for (int i=0; i<N; i++)
	{
		voxels_m[i] = new double*[N];
		voxels_f[i] = new double*[N];
		for (int j=0; j<N; j++)
		{
			voxels_m[i][j] = new double[H];
			voxels_f[i][j] = new double[H];
		}
	}


	time(&t1);
	for (Time=1; true; Time++)
	{
		if (Time/300*300 == Time)
		{
			Excite(Vm);
			Excite(Vm);
			//Excite(Vm);
		}
		

	//	Get_f_sum(Vm,Vf,fs_m,fs_f,fs_e,f_sum,L_Vm,L_Vf);// - this is the serial version of the function calculaing
		// the total currents entering the cells. Now it is substituted by the parallel_for:GetTotalCurrent
		parallel_for(blocked_range<int>(0,N*N*H),GetTotalCurrent(Vm,Vf,fs_m,fs_f,fs_e,f_sum,L_Vm,L_Vf));
		
		
		//SolvePoisson(f_sum,N,H,h,sm_y+sf_y+se_y,sm_x+sf_x+se_x,sm_z+sf_z+se_z,buff,Fe);// - this is the serial poisson equation solver.
		// the do-while cycle and two parallel_for sections inside produce the parallel solution of the 3 dimensional poisson equation
		//using iterative scheme with the accuracy 0.01;
		do
		{
			parallel_for(blocked_range<int>(0,N*N*H),PoissonSolver(f_sum,sm_y+sf_y+se_y,sm_x+sf_x+se_x,sm_z+sf_z+se_z,buff,Fe));
			parallel_for(blocked_range<int>(0,N*N*H),SolvePoissonUpdate(buff,Fe));
		}
		while (error > 0.01);
		
		//Get_external_currents(Vm,Vf,Fe,fs_m,fs_f,Im_e,If_e,L_Vm,L_Vf,L_Fe1,L_Fe2);// - this is the serial function calculating the ionic current
		//going through the each cell. The parallel version of it is the parallel_for:GetExternalCurrents.
		parallel_for(blocked_range<int>(0,N*N*H),GetExternalCurrents(Vm,Vf,Fe,fs_m,fs_f,Im_e,If_e,L_Vm,L_Vf,L_Fe1,L_Fe2));
	
		
		//ANNOTATE_SITE_BEGIN(MySite1);
		//for (i=0; i<N*N*H; i++)
		//{
		//	ANNOTATE_TASK_BEGIN(Task1);
		//	OdeSolve_myocyte(Vm[i],mG[i],hG[i],jG[i],dG[i],fG[i],XG[i],Cai[i]);		
		//	OdeSolve_fib(i,Vf,FB);
		//	ANNOTATE_TASK_END(Task1);
		//}
		//ANNOTATE_SITE_END(MySite1);
		//ANNOTATE_SITE_BEGIN(MySite2);
		//for (i=0; i<N*N*H; i++)           //this commented section integrates the cells' individuals dynamics over the time step "dt"
		//{
		//	ANNOTATE_TASK_BEGIN(Task2);//(first "for" cycle). Then the coupling between the cells is added based on the previously calculated
		//	Vm[i]+=dt*Im_e[i];			  //currents (the second "for" cycle).
		//	Vf[i]+=dt*If_e[i];			  //The same is done by the parallel_for:Psolve
		//	ANNOTATE_TASK_END(Task2);
		//}
		//ANNOTATE_SITE_END(MySite2);
		parallel_for(blocked_range<int>(0,N*N*H),Psolve(Vm,mG,hG,jG,dG,fG,XG,Cai,Vf,FB,Im_e,If_e,dVm,dVf));

		

		if ((Time/Redisplay_time*Redisplay_time == Time) &&(Time >0))
		{
			DrawEvent = true;
			//printf("time=%i\n",Time);
		}
		//else 
			

		if (DrawEvent&&DrawReady)
		{
			DrawEvent = false;
			//vertices_m.clear();
			//getVoxels_d(Vm,N, N, H, voxels);
			parallel_for(blocked_range3d<int>(0,N,0,N,0,H),getVoxels(Vm,voxels_m));
			parallel_for(blocked_range3d<int>(0,N,0,N,0,H),getVoxels(Vf,voxels_f));

			//vertices_m = runMarchingCubes_d(voxels,dVm, -60.0);
			//vertices_f.clear();
			//getVoxels_d(Vf,N, N, H, voxels);
			//vertices_f = runMarchingCubes_d(voxels,dVf, -60.0);
			//printf("vertices_m.size=%i\t",vertices_m.size());
			//printf("vertices_f.size=%i\t",vertices_f.size());
			DrawReady = false;
			glutPostRedisplay();
			
		}

		time(&t2);
		speed = double(Time)/double(t2-t1);
	}

	
	return 1.;
}

void OdeSolve_myocyte(double &Vm, double &mG, double &hG, double &jG, double &dG, double &fG, double &XG, double &Cai, double &dVm)
{
	//this function performs integration of the myocite cell over the time step "dt" possibly subdividing this time step inte several substeps
	double vd;


	dVm = Vm;


		vd=VFunction(Vm,mG,hG,jG,dG,fG,XG,Cai);
		Cai += dt*CaiFunction(Cai,dG,fG,Vm);
		mG   = mFunction(Vm,mG,dt);
		hG   = hFunction(Vm,hG,dt);
		jG   = jFunction(Vm,jG,dt);
		dG   = dFunction(Vm,dG,dt);
		fG   = fFunction(Vm,fG,dt);
		XG   = XFunction(Vm,XG,dt);
		Vm  += dt*vd;

	dVm = (Vm-dVm)/dt;
}

inline int Substeps(double &vd)
{
	// subdivides the time step "dt" into "k" substeps proportionally to the value of the first time derivative of the cell voltage "vd"
	const int kmax=100;
	int k;
	const int k0=vd>0. ? 5 : 1;
 	k=k0+(int)fabs(vd);
	return k<kmax ? k : kmax;
}



void OdeSolve_fib(int i, double *V, Fibroblast *FB)
{
	//integrates the fibroblast dynamics over the time step "dt"
	double dV, dC0, dC1, dC2, dC3, dC4, dO;	
	dV = dt*Vf_function(V[i],FB[i].O_shkr);
	dC0 = dt*C0_function(FB[i].C0,FB[i].C1,V[i]);
	dC1 = dt*C1_function(FB[i].C0,FB[i].C1,FB[i].C2,V[i]);
	dC2 = dt*C2_function(FB[i].C1,FB[i].C2,FB[i].C3,V[i]);
	dC3 = dt*C3_function(FB[i].C2,FB[i].C3,FB[i].C4,V[i]);
	dC4 = dt*C4_function(FB[i].C3,FB[i].C4,FB[i].O_shkr,V[i]);
	dO = dt*O_function(FB[i].C4,FB[i].O_shkr);

	V[i] += dV;
	FB[i].C0 += dC0;
	FB[i].C1 += dC1;
	FB[i].C2 += dC2;
	FB[i].C3 += dC3;
	FB[i].C4 += dC4;
	FB[i].O_shkr += dO;
}

void Get_f_sum(double *Vm, double *Vf, double *fs_m, double *fs_f, double *fs_e, double *f_sum, double *L_Vm, double *L_Vf)
{
	// serial function for calculation the total interdomain current
	// NOT USED IN THE PARALLEL CASE
	int i;
	int n1, n2, n3;
	int rn, ln, un, dn, tn, bn;
	for (i=0; i<N*N*H; i++) //inner grid points
	{
			n3 = i/(N*N);
			n1 = (i-N*N*n3)/N;
			n2 = i-n1*N-N*N*n3;
			
			if (n1 > 0) ln = n3*N*N + (n1-1)*N + n2;
			else ln = i;

			if (n1 < N-1) rn = n3*N*N + (n1+1)*N + n2;
			else rn = i;

			if (n2 > 0) un = n3*N*N + n1*N + n2 - 1;
			else un = i;
			
			if (n2 < N-1) dn = n3*N*N + n1*N + n2 + 1;
			else
				dn = i;

			if (n3 > 0)
				tn = (n3-1)*N*N + n1*N + n2;
			else 	
				tn = i;

			if (n3 < H-1) bn = (n3+1)*N*N + n1*N + n2;
			else bn = i;
			L_Vm[i] = 1e-3*(sm_x*(Vm[rn]+Vm[ln]-2.*Vm[i])+sm_y*(Vm[un]+Vm[dn]-2.*Vm[i])+sm_z*(Vm[tn]+Vm[bn]-2.*Vm[i]))/(h*h);
			L_Vf[i] = 1e-3*(sf_x*(Vf[rn]+Vf[ln]-2.*Vf[i])+sf_y*(Vf[un]+Vf[dn]-2.*Vf[i])+sf_z*(Vf[tn]+Vf[bn]-2.*Vf[i]))/(h*h);
	}
	for (i=0; i<N*N*H; i++)
		f_sum[i] = -L_Vm[i]-L_Vf[i]-fs_m[i]-fs_f[i]-fs_e[i];
}


void Get_external_currents(double *Vm, double *Vf, double *Fe, double *fs_m, double *fs_f, double *Im_e, double *If_e, double *L_Vm, double *L_Vf, double *L_Fe1, double *L_Fe2)
{
	//Calculates the external currents for fibroblasts and myocites
	// NOT USED IN THE PARALLEL CASE
	int i;
	int n1, n2, n3;
	int rn, ln, un, dn, vn, bn, s;
	
	for (i=0;i<N*N*H; i++)
	{
		n3 = i/(N*N);
		n1 = (i-N*N*n3)/N;
		n2 = i-n1*N-N*N*n3;

		n3++;n1++;n2++;

		rn = n3*(N+2)*(N+2) + (n1+1)*(N+2) + n2;
		ln = n3*(N+2)*(N+2) + (n1-1)*(N+2) + n2;
		un = n3*(N+2)*(N+2) + n1*(N+2) + n2 - 1;
		dn = n3*(N+2)*(N+2) + n1*(N+2) + n2 + 1;
		vn = (n3-1)*(N+2)*(N+2) + n1*(N+2) + n2;
		bn = (n3+1)*(N+2)*(N+2) + n1*(N+2) + n2;
		s = n3*(N+2)*(N+2) + n1*(N+2) + n2;

		L_Fe1[i] = 1e-3*(sm_x*(Fe[rn]+Fe[ln]-2.*Fe[s])+sm_y*(Fe[un]+Fe[dn]-2.*Fe[s])+sm_z*(Fe[vn]+Fe[bn]-2.*Fe[s]))/(h*h);
		L_Fe2[i] = 1e-3*(sf_x*(Fe[rn]+Fe[ln]-2.*Fe[s])+sf_y*(Fe[un]+Fe[dn]-2.*Fe[s])+sf_z*(Fe[vn]+Fe[bn]-2.*Fe[s]))/(h*h);
		Im_e[i] = (L_Vm[i]+L_Fe1[i]+fs_m[i]-Betta_myo_fib*(Vm[i]-Vf[i])/(1000.*R_myo_fib))/Betta_myo; 
		If_e[i] = (L_Vf[i]+L_Fe2[i]+fs_f[i]+Betta_myo_fib*(Vm[i]-Vf[i])/(1000.*R_myo_fib))/Betta_fib; 
	}
}


void SolvePoisson(double *f_sum, int N, int H, double h, double ssx, double ssy, double ssz, double *buff, double *Fe)
{
	//performs the serial sulution of the 3dimensional poisson equation using iterative schem with error level of 0.01
	double error, er;
	int i;
	int n1,n2,n3;
	int rn, ln, un, dn, vn, bn;
	//ANNOTATE_SITE_BEGIN(MySite2);
	do
	{
		error = 0;
		
		for (i=0;i<N*N*H; i++)
		{
		//	ANNOTATE_TASK_BEGIN(taskP);
			n3 = i/(N*N);
			n1 = (i-N*N*n3)/N;
			n2 = i-n1*N-N*N*n3;

			n3++;n1++;n2++;

			rn = n3*(N+2)*(N+2) + (n1+1)*(N+2) + n2;
			ln = n3*(N+2)*(N+2) + (n1-1)*(N+2) + n2;
			un = n3*(N+2)*(N+2) + n1*(N+2) + n2 - 1;
			dn = n3*(N+2)*(N+2) + n1*(N+2) + n2 + 1;
			vn = (n3-1)*(N+2)*(N+2) + n1*(N+2) + n2;
			bn = (n3+1)*(N+2)*(N+2) + n1*(N+2) + n2;
			buff[i] = (ssz*(Fe[vn]+Fe[bn])+ssx*(Fe[rn]+Fe[ln])+ssy*(Fe[un]+Fe[dn])-f_sum[i]*h*h)/(2*(ssx+ssy+ssz));
			
			er = fabs(buff[i]-Fe[n3*(N+2)*(N+2)+n1*(N+2)+n2]);
			if (er > error) error = er;
			//ANNOTATE_TASK_END(taskP);
		}

		for (i=0;i<N*N*H; i++)
		{
		//	ANNOTATE_TASK_BEGIN(taskP1);
			n3 = i/(N*N);
			n1 = (i-N*N*n3)/N;
			n2 = i-n1*N-N*N*n3;
			n3++;n1++;n2++;
			Fe[n3*(N+2)*(N+2)+n1*(N+2)+n2] = buff[i];
		//	ANNOTATE_TASK_END(taskP1);
		}
	}
	while (error > 0.01);
	//ANNOTATE_SITE_END(MySite2);
}

