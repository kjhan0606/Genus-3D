#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>




int main(int argc, char **argv){
	int nx,ny,nz,nsub,cubicsize;
	int cubesizex,cubesizey,cubesizez;
	int nsubx,nsuby,nsubz;



	nx = ny = nz = atoi(argv[1]);
	nsub = atoi(argv[2]);

	cubicsize = (nx/nsub)*(nx/nsub)*(nx/nsub); 
	cubesizez = cubesizey = pow(cubicsize/2.,0.333333333333L); 
	cubesizex = (double) cubicsize/(double)cubesizey/(double)cubesizez; 
	nsubx = nx/cubesizex; 
	nsuby = ny/cubesizey; 
	nsubz = nz/cubesizez;
	
	double vol = cubesizex*cubesizey*cubesizez;
	vol = vol * pow(10815.L/(double) nx,3.L);
	printf("subplate volume = %g\n",vol);



}
