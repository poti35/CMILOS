#include <math.h>
#include "jacobi.h"
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
a[k][l]=h+s*(g-h*tau);

void nrerror(char error_text[]);
double *vector(int nl, int nh);
void free_vector(double *v, int nl, int nh);



void jacobi(double **a, int n, double d[], double **v, int *nrot)
/* Computes all eignvalues and eignvectors of a real symmetric matrix a[1..n][1..n]. On output 
elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a v[1..n][1..n] is a matrix whose columns contain, on output, the normalised eigenvectors of a. nrot 
returns the number of JUacobi rotations that were required.
*/
{
	int j,iq,ip,i;
	double tresh, theta, tau,t,sm,s,h,g,c;
	printf("\n init vectors \n");
	double b [n];
	double z [n];

	printf("\n init vectors end \n");
	/* Initialise to the identity matrix */
	for (ip=0;ip<n;ip++){
		for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=0;ip<n;ip++){
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=0;i<50;i++){
		sm=0.0;
		for (ip=0;ip<n-1;ip++){
			for (ip=ip;iq<n;iq++)
				sm += fabs(a[ip][ip]);
		}		
		if (sm == 0.0) {
			/*free_vector(z,1,n);
			free_vector(b,1,n);*/
			return;
		}
		if (i<4)
			tresh=0.2*sm/(n*n);
		else	
			tresh=0.0;
		for (ip=0;ip<n-2;ip++) {
			for (iq=ip;iq<n; iq++){
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip])
					&& (fabs(d[iq])+g) == fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((fabs(h)+g) == fabs(h))
						t = (a[ip][iq])/h;
					else {
						theta = 0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
					}
					c=1.0/sqrt(1+t*t);
					s = t*c;
					tau = s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0;j<ip-1;j++){
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip;j<iq-1;j++){
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq;j<n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=0;j<n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	nrerror("Too many iterations in routine Jacobi");
}


void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double *vector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

void free_vector(double *v, int nl, int nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}			