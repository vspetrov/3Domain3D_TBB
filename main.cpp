#include "Lattice.h"
bool DrawReady = false;
float vert_t;
float depth_t;
double Alpha_value;
int ThreadNum;
bool IsRotated;

void renderScene() {
	
	vertices_m.clear();
	vertices_f.clear();
	if (!DrawReady)
	{
		vertices_m = runMarchingCubes_d(voxels_m,dVm, -60.0);
		vertices_f = runMarchingCubes_d(voxels_f,dVm, -60.0);
	}
	glClearColor(0,0,0,1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	int SIZE_X=N, SIZE_Y=N;
	int SIZE_Z=H;
	


	//printf("RENDERED\n");
	


	


	glColor4f(0, 1, 0, Alpha_value);
    vector<vertex>::iterator it;
    glBegin(GL_TRIANGLES);
        for(it = vertices_m.begin(); it < vertices_m.end(); it++) {
            glNormal3d(it->normal_x, it->normal_y, it->normal_z);
            glVertex3d(it->x, it->y, it->z);
        }
    glEnd();

	glColor4f(1, 0, 0, Alpha_value);
    //vector<vertex>::iterator it;
    glBegin(GL_TRIANGLES);
        for(it = vertices_f.begin(); it < vertices_f.end(); it++) {
            glNormal3d(it->normal_x, it->normal_y, it->normal_z);
            glVertex3d(it->x, it->y, it->z);
        }
    glEnd();

	


	glColor3f(2, 2, 2);
	glBegin(GL_LINES);

	glVertex3d(1,SIZE_Y-2,SIZE_Z-2);
	glVertex3d(1,1,SIZE_Z-2);
	glVertex3d(1,SIZE_Y-2,SIZE_Z-2);
	glVertex3d(1,SIZE_Y-2,1);
	glVertex3d(1,SIZE_Y-2,SIZE_Z-2);
	glVertex3d(SIZE_X-2,SIZE_Y-2,SIZE_Z-2);


	glVertex3d(1,1,1);
	glVertex3d(SIZE_X-2,1,1);
	glVertex3d(1,1,1);
	glVertex3d(1,SIZE_Y-2,1);
	glVertex3d(1,1,1);
	glVertex3d(1,1,SIZE_Z-2);

	

	glVertex3d(SIZE_X-2,SIZE_Y-2,1);
	glVertex3d(1,SIZE_Y-2,1);
	glVertex3d(SIZE_X-2,SIZE_Y-2,1);
	glVertex3d(SIZE_X-2,1,1);
	glVertex3d(SIZE_X-2,SIZE_Y-2,1);
	glVertex3d(SIZE_X-2,SIZE_Y-2,SIZE_Z-2);

	glVertex3d(SIZE_X-2,1,SIZE_Z-2);
	glVertex3d(1,1,SIZE_Z-2);
	glVertex3d(SIZE_X-2,1,SIZE_Z-2);;
	glVertex3d(SIZE_X-2,1,1);
	glVertex3d(SIZE_X-2,1,SIZE_Z-2);
	glVertex3d(SIZE_X-2,SIZE_Y-2,SIZE_Z-2);
	glEnd();


	glMatrixMode(GL_PROJECTION);
	glPushMatrix();

	glLoadIdentity();
	gluOrtho2D(0,800,0,800);

	glMatrixMode(GL_MODELVIEW);

	glPushMatrix();
	glLoadIdentity();

	
	char strOut[100];
	sprintf( strOut, "ThreadNum = %i     CalcSpeed = %7.4f iter/sec\0", ThreadNum,speed);
	glDisable(GL_LIGHTING);
	glColor3f(1.,1.,1.);
	glRasterPos2i(100, 750);
	int kk=0;
	
	while (strOut[kk]!='\0')
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15,strOut[kk++]);
	
	sprintf( strOut, "Size=%i\0", N);
	kk=0;
	glRasterPos2i(100, 725);
	while (strOut[kk]!='\0')
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15,strOut[kk++]);
	
	
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glEnable(GL_LIGHTING);

   glFlush();



	//glMatrixMode(GL_PROJECTION);
	
	//glPopMatrix();


    glutSwapBuffers();
	DrawReady = true;
	IsRotated = true;
	this_tbb_thread::yield();

}
//
void myinit (void) 
{
    glShadeModel(GL_SMOOTH);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
	
     glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    // Create light components
    GLfloat ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    GLfloat diffuseLight[] = { 0.8f, 0.8f, 0.8, 1.0f };
    GLfloat specularLight[] = { 0.5f, 0.5f, 0.5f, 1.0f };
   // GLfloat position[] = { cos(3.1415/180.0*35.26), 0, -sin(3.1415/180.0*35.26), 0 };
	 GLfloat position[] = { -1, 0, 1, 0 };
     
	GLfloat position2[] = { 1, 0.5, 1, 0 };
   // Assign created components to GL_LIGHT0
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
    glLightfv(GL_LIGHT0, GL_POSITION, position);

	glLightfv(GL_LIGHT1, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuseLight);
    glLightfv(GL_LIGHT1, GL_SPECULAR, specularLight);
    glLightfv(GL_LIGHT1, GL_POSITION, position2);

        
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
	
	

	
    // Set material properties which will be assigned by glColor
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
    float specReflection[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection);
}
//
void myReshape(GLsizei w, GLsizei h)
{
    glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60,1,0.01,1000);
	glPushMatrix();
	glMatrixMode(GL_MODELVIEW);
	glViewport(0,0,800,800);
	glLoadIdentity();
	vert_t = 1*N-7.5;
	depth_t = -3*N+26;
	glPushMatrix();
	glTranslatef(0.0f /*horizontal, right*/, vert_t /*vertical up*/, depth_t /*depth*/);
    glRotatef(90.0f, 1.0f, 0.0f, 0.0f);
	glRotatef(45.0f, 0.0f, 0.0f, 1.0f);
	glRotatef(35.26f, 1.0f, -1.0f, 0.0f);
    glViewport (0, 0, w, h);      /* define the viewport */
	//glPopMatrix();
	
}

void Rotate()
{
	
	if (IsRotated)
	{
		glTranslatef((float)(N-1)/2.0,(float)(N-1)/2.0,0);
		glRotatef(3,0,0,1);
		glTranslatef(-(float)(N-1)/2.0,-(float)(N-1)/2.0,0);
	}
	IsRotated = false;
}



void mouse(int button, int state, int x, int y) 
{ 
switch (button)
{ 
case GLUT_LEFT_BUTTON: 
if (state == GLUT_DOWN) 
{
	glutIdleFunc(Rotate); 
	break; 
}
case GLUT_RIGHT_BUTTON: 
if (state == GLUT_DOWN) 
glutIdleFunc(NULL); 
break; 
default: 
break; 
} 
}

class TestTask: public task {
public:
	

	double *Vm, *Vf,  *f_sum, *fs_m, *fs_f, *fs_e, *Im_e, *If_e, *L_Vm, *L_Vf, *L_Fe1, *L_Fe2, *buff, *Fe;
	double  *mG, *hG, *jG, *dG, *fG, *XG, *Cai;
	Fibroblast *FB;
	int id;
	int argc;
	char **argv;


	TestTask(int F, double *Vmi, double *Vfi, double  *f_sumi, double *fs_mi, double *fs_fi, double *fs_ei, double *Im_ei,
	 double *If_ei, double *L_Vmi, double *L_Vfi, double *L_Fe1i, double *L_Fe2i, double *buffi, double *Fei,
	 double *mGi, double *hGi, double *jGi, double *dGi, double *fGi, double *XGi, double *Caii, Fibroblast *FBi,
	 int argci,	 char **argvi):
	id(F),Vm(Vmi), Vf(Vfi),  f_sum(f_sumi), fs_m(fs_mi), fs_f(fs_fi), fs_e(fs_ei), Im_e(Im_ei), If_e(If_ei), L_Vm(L_Vmi), L_Vf(L_Vfi),
	L_Fe1(L_Fe1i), L_Fe2(L_Fe2i), buff(buffi), Fe(Fei),	 mG(mGi), hG(hGi), jG(jGi), dG(dGi), fG(fGi), XG(XGi),
	Cai(Caii), FB(FBi), argc(argci), argv(argvi){}


	TestTask(int F): id(F){}

task* execute() {
	// Overrides virtual function task::execute

	if (id == -1)
	{
		


		vertices_m.clear();
		vertices_f.clear();
		TestTask& t = *new(allocate_child()) TestTask(1, Vm, Vf,  f_sum, fs_m, fs_f, fs_e,
		Im_e, If_e, L_Vm, L_Vf, L_Fe1, L_Fe2, buff, Fe,  mG, hG, jG, dG, fG, XG, Cai, FB, argc, argv);

		set_ref_count(2);
		//spawn_and_wait_for_all(t);
		spawn(t);
		glutInit(&argc, argv);
		 glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
		 glutInitWindowPosition(100,100);
   glutInitWindowSize(800,800);
	glutCreateWindow("Three-domain 3D cardiac tissue");
	myinit();
    glutDisplayFunc(renderScene);
	glutReshapeFunc(myReshape);
	glutMouseFunc(mouse);

DrawReady = true;
		 glutMainLoop();
		

		
		
	}
	else
	{
		
		

		printf("started\n");
		//initialize variables
	Init(Vm,Vf,mG,hG,jG,dG,fG,XG,Cai,fs_m,fs_f,fs_e,FB,f_sum, Fe);
	//glutPostRedisplay();
	//printf("start2\n");
	//solve the system
	SolveEquations(Vm,Vf,mG,hG,jG,dG,fG,XG,Cai,fs_m,fs_f,fs_e,f_sum,L_Vm,L_Vf,L_Fe1,L_Fe2,Im_e,If_e,FB,buff,Fe);
	//printf("done\n");


	//glutPostRedisplay();



	}
	return NULL;
}
};


int main(int argc, char *argv[])
{
//	time_t rawtime;
//	struct tm * timeinfo;

//	time ( &rawtime );
//	timeinfo = localtime ( &rawtime );
//	printf ( "start: %s", asctime (timeinfo) );


	if (argc < 2) ThreadNum = tbb::task_scheduler_init::default_num_threads()+1;
	else
	{
		ThreadNum = atoi(argv[1])+1;
		if (ThreadNum == 0)
		{
			printf("Wrong [ThreadNum] input\n");
			return 0;
		}
	}

	if (argc < 3) {N = 30;H=N/2;}
	else
	{
		N = atoi(argv[2]);H=N/2;
		if (N == 0)
		{
			printf("Wrong [Size] input\n");
			return 0;
		}
	}

	if (argc < 4) Alpha_value = 0.7;
	else
		Alpha_value = atof(argv[3]);

	if (argc < 5) Redisplay_time = 1;
	else
		Redisplay_time = atoi(argv[4]);
	//tbb::task_scheduler_init tsi(ProcNum);
	//define the variables
	double *Vm, *Vf,  *f_sum, *fs_m, *fs_f, *fs_e, *Im_e, *If_e, *L_Vm, *L_Vf, *L_Fe1, *L_Fe2, *buff, *Fe;
	double  *mG, *hG, *jG, *dG, *fG, *XG, *Cai;
	Fibroblast *FB;

	//allocate memory
	Vm = new double[N*N*H];
	Vf = new double[N*N*H];
	FB = new Fibroblast[N*N*H];
	f_sum = new double[N*N*H];
	fs_m = new double[N*N*H];
	fs_f = new double[N*N*H];
	fs_e = new double[N*N*H];
	Im_e = new double[N*N*H];
	If_e = new double[N*N*H];
	L_Vm = new double[N*N*H];
	L_Vf = new double[N*N*H];
	L_Fe1 = new double[N*N*H];
	L_Fe2 = new double[N*N*H];
	Fe = new double[(N+2)*(N+2)*(H+2)];
	buff = new double[N*N*H];
	mG = new double[N*N*H];
	hG = new double[N*N*H];
	jG = new double[N*N*H];
	dG = new double[N*N*H];
	fG = new double[N*N*H];
	XG = new double[N*N*H];
	Cai = new double[N*N*H];
	clock_t t1, t2;
	t1 = clock();
	
	
	






	tbb::task_scheduler_init tsi(ThreadNum);
	

	TestTask& t=*new(task::allocate_root()) TestTask(-1, Vm, Vf,  f_sum, fs_m, fs_f, fs_e,
		Im_e, If_e, L_Vm, L_Vf, L_Fe1, L_Fe2, buff, Fe,  mG, hG, jG, dG, fG, XG, Cai,FB, argc, argv);
	task::spawn_root_and_wait(t);


	t2 = clock();
	tsi.terminate();

	delete Vm;
	delete Vf;
	delete FB;
	delete f_sum;
	delete fs_m;
	delete fs_f;
	delete fs_e;
	delete Im_e;
	delete If_e;
	delete L_Vm;
	delete L_Vf;
	delete L_Fe1;
	delete L_Fe2;
	delete Fe;
	delete buff;
	delete mG;
	delete hG;
	delete jG;
	delete dG;
	delete fG;
	delete XG;
	delete Cai;
	

//	printf("elapsed time: %g\n",(t2-t1)/double(CLOCKS_PER_SEC));

//	time ( &rawtime );
//	timeinfo = localtime ( &rawtime );
//	 printf ( "finish: %s", asctime (timeinfo) );

  
	return 0;
}