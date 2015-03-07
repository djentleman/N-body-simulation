#include <stdlib.h>
#include <stdio.h>
#include <windows.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>
#include <time.h>




static GLfloat theta[]={0.0,0.0,0.0};
static GLdouble viewer[]={0.0,0.0,5.0};

static double G = 0.000001;//0.00000000006674; // newtons gravitational constant

static int N = 100;
static int DELAY = 0; // delay between frames (microseconds)
static double dT = 0.02; // change it rate of time

static int p = 2; // number of dimensions, either 2 or 3

int iter = 0;

int miter = 500; // max iterations
bool bound = true; // if bound, stop when iter == miter, else run forever

static int height = 800;
static int width = 800;

double M = 1000000; // maximum mass, only works if rM
double GC = 10000000; // centeral galactic mass
double R = 6; // body spawn radius

bool rM = false; // random masses
bool rV = false; // random initial velocity
bool ex = false; // expansion
bool rot = false; // rotation
bool gc = false; // galactic center
bool clu = false; // collisio2

clock_t start;
clock_t end; // for measuring seconds elapsed

bool img = false;

// collision speeds
double cY = -0.0;
double cX = -10;

double IR = 0.9; // imbalance of clusters

double er = 1; // expansion rate
double rr = 1; // rotation rate

int T = 8; // thread count

bool gauss = true; // generate randoms from a gaussian distribution

bool follow = true; // camera follow average point
bool lighting = false;

bool cmd = true; // allows command line input from user



// position arrays
double* px;
double* py;
double* pz;
// velocity arrays
double* vx;
double* vy;
double* vz;
// mass array
double* m;




double getMean(int n, double* arr)
{
	int i;
	double sum = 0;
	for (i = 0; i < n; i++)
	{
		sum += arr[i];
	}
	return sum/n;
}

double getMinkowskiDistance(int i, int j)
{
	double norm = 0;
	norm += (double)pow((double)fabs(px[i] - px[j]), (double)p);
	//printf("%f %f %f, %f %f %f -> %e\n", px[i], py[i], pz[i], px[j], py[j], pz[j], norm);
	norm += pow(fabs(py[i] - py[j]), p);
	if (p == 3)
		norm += pow(fabs(pz[i] - pz[j]), (double)p);
	
	norm = pow(norm, (double)1/(double)p);
	//printf("%f %f %f, %f %f %f -> %e\n", px[i], py[i], pz[i], px[j], py[j], pz[j], norm);
	return pow(norm, p);
}

double getForce(int i, int j, int d)
{
	// gets effect j has on i
	// F = G*M*M (i - j) / getMinkowskiDistance(i, j)
	
	if (d == 0)
	{
		// x
		return G*m[i]*m[j]*(px[i] - px[j]) / getMinkowskiDistance(i, j);
	}
	else if (d == 1)
	{
		// y
		return G*m[i]*m[j]*(py[i] - py[j]) / getMinkowskiDistance(i, j);
	}
	else
	{
		// z;
		return G*m[i]*m[j]*(pz[i] - pz[j]) / getMinkowskiDistance(i, j);
	}



}

void updateVelocity(int i)
{
	
	int j;
	double Fx = 0;
	double Fy = 0;
	double Fz = 0;
	for (j = 0; j < N; j++)
	{
		if (i != j)
		{
			Fx += getForce(i, j, 0);
			Fy += getForce(i, j, 1);
			if (p == 3)
				Fz += getForce(i, j, 2);
			
			//printf("%e, %e, %e\n", Fx, Fy, Fz);
		}
	}
	// update tx, ty, tz
	// F = MA -> A = F/M
	// A = dV -> v(n+1) = v(n) + A, simple :)
	double Ax = Fx/m[i];
	double Ay = Fy/m[i];
	double Az = 0;
	if (p == 3)
		Az = Fz/m[i];
	
	vx[i] = vx[i] - (Ax * dT);
	vy[i] = vy[i] - (Ay * dT);
	if (p == 3) // else will stay 0
		vz[i] = vz[i] - (Az * dT);
}

void init()
{
    // set up open MP
	omp_set_dynamic(0);
	omp_set_num_threads(T);
	
	
	px = malloc(sizeof(double) * N);
	py = malloc(sizeof(double) * N);
	pz = malloc(sizeof(double) * N);
	
	vx = malloc(sizeof(double) * N);
	vy = malloc(sizeof(double) * N);
	vz = malloc(sizeof(double) * N);
	
	
	m = malloc(sizeof(double) * N);
	
	srand(time(0));
	int i;
	double rnd;
	for (i = 0; i < N; i++)
	{
		rnd = ((double)rand()/RAND_MAX);
		if (clu)
		{
			if (rnd > IR)
			{
				// cluster 1
				px[i] = (((double)rand()/RAND_MAX)-0.5)*2;
				py[i] = (((double)rand()/RAND_MAX)-0.5)*(sqrt(1-pow(px[i], 2)))*2;
			}
			else
			{
				// cluster 2
				px[i] = (((double)rand()/RAND_MAX)-0.5+R);
				py[i] = (((double)rand()/RAND_MAX)-0.5)*((sqrt(pow(0.5, 2) - pow(px[i]-R, 2)))*2)+R;
			}
			
			if (p == 2)
			{
				pz[i] = 0;
			}
			else
			{
				pz[i] = (((double)rand()/RAND_MAX)-0.5);
			}	
		}
		else
		{
			
			
			
			px[i] = (((double)rand()/RAND_MAX)-0.5)*R;
			// get max possible value of y
			//a^2+b^2 = c^2
			//b = sqrt((R/2)^2-px[i]^2)
			py[i] = (((double)rand()/RAND_MAX)-0.5)*(sqrt(pow(R/2, 2)-pow(px[i], 2)))*2; // circular 
			
			
			if (p == 2)
			{
				pz[i] = 0;
			}
			else
			{
				pz[i] = (((double)rand()/RAND_MAX)-0.5)*R;
			}	
			
			if (rot && i == 0 && gc)
			{
				px[0] = 0;
				py[0] = 0;
				pz[0] = 0;
			}
			
			
			
		}
		
		
		if (rV)
		{
			vx[i] =(((double)rand()/RAND_MAX)-0.5);
			vy[i] =(((double)rand()/RAND_MAX)-0.5);
			if (p == 2 || !rV)
			{
				vz[i] = 0;
			}
			else
			{
				vz[i] = (((double)rand()/RAND_MAX)-0.5);
			}
		}
		else if (clu)
		{
			if (rnd > IR)
			{
				vx[i] = 0;
				vy[i] = 0;
			}
			else
			{
				vx[i] = -cX;
				vy[i] = -cY;
			}
		}
		else if (ex)
		{
			vx[i] = px[i]*er;
			vy[i] = py[i]*er;
			vz[i] = pz[i]*er;
		}
		else if (rot)
		{
			vx[i] = py[i]*rr;
			vy[i] = -px[i]*rr;
			vz[i] = 0;
		}
		else
		{
			vx[i] = 0;
			vy[i] = 0;
			vz[0] = 0;
		}
		
		
		//printf("Position: %f %f %f", px[i], py[i], pz[i]);
		//printf("	Velocity: %f %f %f\n", vx[i], vy[i], vz[i]);
		if (rot && i == 0 && gc)
		{
			m[0] = GC;
		}
		else if (rM)
		{
			m[i] = (((double)rand()/RAND_MAX))*M;
		}
		else
		{
			m[i] = M/10;
		}
	}
}

void plot()
{
	// volume of sphere: v = 4pir^3/3
	// cuberoot(3v/4pi) = r
	
	int i;
	for (i = 0; i < N; i++)
	{
		glTranslatef(px[i], py[i], pz[i]);
		glutSolidSphere( cbrt(m[i])/1000, 20, 16);
		glTranslatef(-px[i], -py[i], -pz[i]);
	}
}

void update()
{
	// update velocities based on other points
	int i;
	
	#pragma omp parallel
	{
		#pragma omp for
		for (i = 0; i < N; i++)
		{
			updateVelocity(i);
		}
	}
	
	
	#pragma omp parallel
	{
		#pragma omp for
		for (i = 0; i < N; i++)
		{
			px[i] += vx[i] * dT;
			py[i] += vy[i] * dT;
			pz[i] += vz[i] * dT;
		}
	}
	
	
}

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	if (follow)
	{
		gluLookAt(getMean(N, px)+theta[0], getMean(N, py)+theta[1],getMean(N, pz)+10+theta[2]
		,getMean(N, px)+theta[0], getMean(N, py)+theta[1],getMean(N, pz),
		1.0,1.0,1.0);
	}
	else
	{
		gluLookAt(viewer[0]+theta[0], viewer[1]+theta[1],viewer[2]+theta[2],0.0,0.0,0.0, 
		1.0,1.0,1.0);
	}
	//printf("%f %f %f",getMean(N, px), getMean(N, py),getMean(N, pz));
	//printf("%f %f %f",theta[0], theta[1], theta[2]);
	
	if (lighting)
	{
		GLfloat light_position[] = { 3.0, 0.0, 0.0, 1.0 };


		glLightfv(GL_LIGHT0, GL_POSITION, light_position);

		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_DEPTH_TEST);
		
	}
	if (img) {
		
		} // do nothing
	
	printf("Iter: %i\n", ++iter);
	plot();
	glFlush();
	
	glutSwapBuffers();
	update();
	usleep(DELAY);
	if (bound)
	{
		if (iter >= miter)
		{
			end = clock();
			printf("Elapsed: %f\n", ((float)end - (float)start));
			return;
		}
	}
	glutPostRedisplay();
}


void keys(unsigned char key, int x, int y)
{
	if(key == 'a') {theta[0]+=0.1;theta[1]-=0.1;}
	if(key == 'd') {theta[0]-=0.1;theta[1]+=0.1;}
	if(key == 'w') {theta[0]-=0.1;theta[1]-=0.1;}
	if(key == 's') {theta[0]+=0.1;theta[1]+=0.1;}
	if(key == 'q') theta[2]-=0.1;
	if(key == 'e') theta[2]+=0.1;
	display();
}

void myReshape(int w, int h)
{
	glViewport(0,0,w,h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if(w<=h)
        glFrustum(-2.0,2.0,-2.0*(GLfloat)h/(GLfloat)w, 
                 2.0*(GLfloat)h/(GLfloat)w, 2.0,20.0);
	else 
        glFrustum(-2.0,2.0,-2.0*(GLfloat)w/(GLfloat)h, 
2.0*(GLfloat)w/(GLfloat)h, 2.0,20.0);
	glMatrixMode(GL_MODELVIEW);
}

char* read()
{
	char *buffer = malloc(sizeof(char) * 100);
    int read;
    read = fgets(buffer, 100, stdin);
	
    if (-1 != read)
        return buffer;
    return NULL;
}

void getInput()
{
	
	char* in;
	int v;
	float f;
	// gets input from user, updates globals
	printf("Nbody.c UI\n");
	printf("Press Enter To Set Parameters To Default\n\n");
	
	
	
	printf("Number Of Bodies: ");
	v = atoi(read());
	if (v != 0){N = v;}
	
	printf("Number Of Dimensions {2, 3}: ");
	v = atoi(read());
	if (v == 2 || v == 3){p = v;}
	
	printf("Number Of Threads: ");
	v = atoi(read());
	if (v != 0){T = v;}
	
	printf("Gravitational Constant: ");
	f = atof(read());
	if (f != 0){G = f;}
	
	printf("dT (Rate Of Change In Time): ");
	f = atof(read());
	if (f != 0){dT = f;}
	
	printf("Spawn Radius: ");
	f = atof(read());
	if (f != 0){R = f;}
	
	printf("Delay Between Frames (micro seconds): ");
	v = atoi(read());
	if (v != 0){DELAY = v;}
	
	printf("Stop After n Iterations? 1=NO, 2=YES: ");
	v = atoi(read());
	if (v == 2)
	{
		bound = true;
		printf("How Many Iterations? ");
		v = atoi(read());
		if (v != 0){miter = v;}
		
	}
	else if (v == 1)
	{
		bound = false;
	}
	
	printf("Show lighting on bodies? 1=NO, 2=YES: ");
	v = atoi(read());
	if (v == 2)
	{
		lighting = true;
	}
	else if (v == 1)
	{
		lighting = false;
	}
	
	printf("Randomly Initialize Masses? 1=NO, 2=YES: ");
	v = atoi(read());
	if (v == 2)
	{
		rM  = true;
	}
	else if (v == 1)
	{
		rM = false;
	}
	
	printf("\nInitial Conditions:\n1=No Initial Velocity\n2=Random Initial Velocity\n3=Rotating Cluster\n4=Expanding Cluster\n5=Cluster Collision\n: ");
	v = atoi(read());
	if (v != 0)
	{
		rV = false; // random initial velocity
		ex = false; // expansion
		rot = false; // rotation
		clu = false; // collisio2
	}
	if (v == 2)
	{
		rV = true; // random initial velocity
	}
	if (v == 3)
	{
		rot = true; // rotation
		// set up rotation params
		
		printf("Angular Velocity: ");
		f = atof(read());
		if (f != 0){rr = f;}
		
		printf("Centeral Galactic Mass? 1=NO, 2=YES: ");
		v = atoi(read());
		if (v == 2)
		{
			gc = true; // galactic center
			printf("Mass Of CGM: ");
			v = atoi(read());
			if (f != 0){rr = v;}
		}
		else if (v == 1)
		{
			GC = false; // galactic center
		}
	}
	if (v == 4)
	{
		ex = true;
		// set up cluster expansion
		printf("Expansion Rate: ");
		f = atof(read());
		if (f != 0){er = f;}
	}
	if (v == 5)
	{
		clu = true;
		// set up cluster collision
		
		printf("Cluster Imbalance: ");
		f = atof(read());
		if (f != 0){IR = f;}
		
		printf("Secondary Cluster X Velocity: ");
		f = atof(read());
		if (f != 0){cX = f;}
		
		printf("Secondary Cluster Y Velocity: ");
		f = atof(read());
		if (f != 0){cY = f;}
		
	}
}



int main(int argc, char **argv)
{
	if (cmd)
	{
		getInput();
	}
	start = clock();
	init();
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
	glutInitWindowSize(height,width);
	char buf[256]; 
	sprintf(buf, "N-Body Simulation, %i Bodies, %i Dimensions, %i Thread(s)", N, p, T);
	glutCreateWindow(buf);
	glutReshapeFunc(myReshape);
	glutDisplayFunc(display);
	glutKeyboardFunc(keys);
	glEnable(GL_DEPTH_TEST);
	glutMainLoop();
	return 0;
}
