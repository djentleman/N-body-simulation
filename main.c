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

static double G = 0.001;// 0.00000000006674; // newtons gravitational constant

static int N = 100;
static int DELAY = 0; // delay between frames (microseconds)
static double dT = 1; // change it rate of time

static int p = 2; // number of dimensions, either 2 or 3

int iter = 0;

double M = 0.1; // maximum mass, only works if rM
double R = 1; // body spawn radius

bool rM = true; // random masses
bool rV = false; // random initial velocity
bool col = true; // collisio2

int T = 8; // thread count

double IR = 0.05; // imbalance of clusters


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

double zoom = 1;

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
		if (col)
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
		}
		else
		{
			px[i] = (((double)rand()/RAND_MAX)-0.5)*R;
			py[i] = (((double)rand()/RAND_MAX)-0.5)*R;
			
			
		}
		if (p == 2)
		{
			pz[i] = 0;
		}
		else
		{
			pz[i] = (((double)rand()/RAND_MAX)-0.5)*R;
		}	
		
		
		if (rV)
		{
			vx[i] =(((double)rand()/RAND_MAX)-0.5)/100;
			vy[i] =(((double)rand()/RAND_MAX)-0.5)/100;
		}
		else
		{
			if (col)
			{
				if (rnd > IR)
				{
					vx[i] = 0;
					vy[i] = 0;
				}
				else
				{
					vx[i] = 0;
					vy[i] = -0.1;
				}
			}
			else
			{
				vx[i] = 0;
				vy[i] = 0;
			}
		}
		
		if (p == 2 || !rV)
		{
			vz[i] = 0;
		}
		else
		{
			vz[i] = (((double)rand()/RAND_MAX)-0.5)/100;
		}
		//printf("Position: %f %f %f", px[i], py[i], pz[i]);
		//printf("	Velocity: %f %f %f\n", vx[i], vy[i], vz[i]);
		
		if (rM)
		{
			m[i] = (((double)rand()/RAND_MAX))*M;
		}
		else
		{
			m[i] = 0.01;
		}
	}
}

void plot()
{
	int i;
	for (i = 0; i < N; i++)
	{
		glTranslatef(px[i], py[i], pz[i]);
		glutSolidSphere(m[i], 20, 16);
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
	gluLookAt(viewer[0], viewer[1],viewer[2],0.0,0.0,0.0, 
1.0,1.0,1.0);
	glRotatef(theta[0],1.0,0.0,0.0);
	glRotatef(theta[1],0.0,1.0,0.0);
	glRotatef(theta[2],0.0,0.0,1.0);
	glScalef(zoom, zoom, zoom);
	
	
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
	if(key == 'x') theta[0]+=1.0;
	if(key == 'y') theta[1]+=1.0;
	if(key == 'z') theta[2]+=1.0;
	if(key == 'a') zoom+=0.1;
	if(key == 'b') zoom-=0.1;
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
	sprintf(buf, "N-Body Simulation, %i Bodies, %i Dimensions", N, p);
	glutCreateWindow(buf);
	glutReshapeFunc(myReshape);
	glutDisplayFunc(display);
	glutKeyboardFunc(keys);
	glEnable(GL_DEPTH_TEST);
	glutMainLoop();
}
