#include <stdlib.h>
#include <windows.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>






static GLfloat theta[]={0.0,0.0,0.0};
static GLdouble viewer[]={0.0,0.0,5.0};

static double G = 0.000000001;//0.00000000006674; // newtons gravitational constant

static int N = 200;
static int DELAY = 0; // delay between frames (microseconds)
static double dT = 0.2; // change it rate of time

static int p = 3; // number of dimensions, either 2 or 3

int iter = 0;

double M = 1000000; // maximum mass, only works if rM
double R = 10; // body spawn radius

bool rM = true; // random masses
bool rV = true; // random initial velocity
bool clu = false; // collisio2

// collision speeds
double cY = 0.8;
double cX = 0.0;

double IR = 0.1; // imbalance of clusters

int T = 8; // thread count

bool follow = true; // camera follow average point
bool lighting = false;



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
	return norm;
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
	
	srand(time(NULL));
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
				px[i] = (((double)rand()/RAND_MAX)-0.5)/2;
				py[i] = (((double)rand()/RAND_MAX)-0.5)/2;
			}
			else
			{
				// cluster 2
				px[i] = (((double)rand()/RAND_MAX)+2);
				py[i] = (((double)rand()/RAND_MAX)+2);
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
			py[i] = (((double)rand()/RAND_MAX)-0.5)*R;
			
			
			if (p == 2)
			{
				pz[i] = 0;
			}
			else
			{
				pz[i] = (((double)rand()/RAND_MAX)-0.5)*R;
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
		else
		{
			vx[i] = 0;
			vy[i] = 0;
		}
		
		
		//printf("Position: %f %f %f", px[i], py[i], pz[i]);
		//printf("	Velocity: %f %f %f\n", vx[i], vy[i], vz[i]);
		
		if (rM)
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
	
	printf("Iter: %i\n", ++iter);
	plot();
	glFlush();
	glutSwapBuffers();
	update();
	usleep(DELAY);
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




int main(int argc, char **argv)
{
	init();
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
	glutInitWindowSize(800,800);
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
