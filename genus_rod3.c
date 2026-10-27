/* 
  mpicc -o genus_rod3  genus_rod3.c -lm -L/home/kjhan/fftw-intel/lib -I/home/kjhan/fftw-intel/include/  -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
 */

#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>

#include "srfftw_mpi.h"
#include "srfftw.h"



#define NSMOOTH 1

float RG[NSMOOTH]={22};
#define NTHRE 1
float rhothre[NTHRE];
int ithre;

int nid,myid;



float rhoc;
float label,vmfrac,delta;


float genus;

double vol,mass;

double mtot,vtot;


long long npsample;
/*
int *indx;
*/


char volmass,pein,dc;

float boxsize;

double aver,std,stdg;


void readnewline(FILE *fp){
	char buff[10000];
	fgets(buff,10000,fp);
}


rfftwnd_mpi_plan plan,iplan;
int local_nz,local_z_start,local_ny_after_transpose;
int local_y_start_after_transpose,total_local_size;


FILE *wp1, *wp2;

float *copyden(float *,int , int , int );
float *Sampling(float *, int ,int ,int ,int ,int ,int );
void summ(int *, int , int , int , int);
void paddinghighlow(int *,int ,int , int );
void CalHighLow(float *, int , int , int , float , int *);
void FILEDUMP(double , double , int , int );
void smoothit(float *,int ,int ,int ,float );
void stat3d(float *, int , int , int , int,int );
void gradx(float *,int , int , int );
float *meshdenread(char *, int *, int *, int *);
void Statistics(float *, int , int , int );

int sorttype( const void *a, const void *b){
	float *aa, *bb;
	aa= (float*)a;
	bb= (float*)b;
	if(*aa < *bb) return -1;
	else if(*aa > *bb) return 1;
	else return 0;
}


void getmeshinfo(char *infile, int *nx, int *ny, int *nz){
	FILE *fp;
	fp = fopen(infile,"r");
	fread(nx,sizeof(int),1,fp);
	fread(ny,sizeof(int),1,fp);
	fread(nz,sizeof(int),1,fp);
	fclose(fp);
}

int nsubx,nsuby,nsubz;
int cubesizex,cubesizey,cubesizez;

int main( int argc, char **argv){
	int i,j,k;
	int ismooth;
	char infile[109];
	int outnum;
	float *den,*iden;
	int nx,ny,nz;
	FILE *fp;


	MPI_Init(&argc, &argv);

	if(argc != 6){
		fprintf(stderr,"./genus indenfile boxsize outnum p/i nsub\n");
		exit(999);
	}


	sprintf(infile,"%s",argv[1]);

	boxsize = atof(argv[2]);
	outnum = atoi(argv[3]);
	pein = *argv[4];
	int nsub = atoi(argv[5]);


	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);



	getmeshinfo(infile,&nx,&ny,&nz);

	int cubicsize = (nx/nsub)*(nx/nsub)*(nx/nsub);
	cubesizez = cubesizey = pow(cubicsize/4.,0.333333333333L);
	cubesizex = (double) cubicsize/(double)cubesizey/(double)cubesizez;
	nsubx = nx/cubesizex;
	nsuby = ny/cubesizey;
	nsubz = nz/cubesizez;

	plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,nz,ny,nx, 
			FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE);
	iplan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,nz,ny,nx, 
			FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE);

	rfftwnd_mpi_local_sizes(plan,&local_nz,&local_z_start, 
			&local_ny_after_transpose,&local_y_start_after_transpose, 
			&total_local_size);

	iden = meshdenread(infile, &nx,&ny, &nz);
	MPI_Barrier(MPI_COMM_WORLD);
	int mx = 2*(nx/2+1);

	double sigma,mean,tsigma,tmean;
	sigma = mean = 0;
	for(k=0;k<local_nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++){
		mean += iden[i+nx*(j+ny*k)];
		sigma += iden[i+nx*(j+ny*k)] * iden[i+nx*(j+ny*k)];
	}
	MPI_Reduce(&mean,&tmean,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&sigma,&tsigma,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(&tmean,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&tsigma,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	tmean = tmean/(double)nx/(double)ny/(double)nz;
	tsigma = tsigma/(double)nx/(double)ny/(double)nz - tmean*tmean;
	tsigma = sqrt(tsigma);
	if(myid==0) printf("mean/sigma = %g %g\n",tmean,tsigma);
	for(i=0;i<NTHRE;i++){
		rhothre[i] = i*tsigma + tmean;
	}

	for(ithre = 0;ithre < NTHRE; ithre++){

	int mmx,mmy,mmz;
	if(myid==0) {
		fp = fopen("contour.inp","r");
		readnewline(fp);
		fscanf(fp,"%d %d %d\n",&mmx,&mmy,&mmz);
		readnewline(fp);
		/*
		fscanf(fp,"%c%c%c  %c%c%c ",&dc,&volmass,&dc,&dc, &pein,&dc);
		*/
		readnewline(fp);
		readnewline(fp);
		readnewline(fp);
	}
	MPI_Bcast(&mmx,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&mmy,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&mmz,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&volmass,1,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&pein,1,MPI_CHAR,0,MPI_COMM_WORLD);



	int sx,sy,sz;
	sx = sy = sz = 300;



	for(ismooth=0;ismooth<NSMOOTH;ismooth++){
		float esmooth = RG[ismooth]*sqrt(2.);
		if(myid==0){
			char out1[190];
			sprintf(out1,"contour%.2d.out%.2d.ithre%.2d",outnum,ismooth,ithre);
			wp2 = fopen(out1,"w");
			sprintf(out1,"plot%.2d.out%.2d.ithre%.2d",outnum,ismooth,ithre);
			wp1 = fopen(out1,"w");
		}


		if(myid==0) printf("before copyden\n");
		den = copyden(iden, nx, ny,nz);
		
		if(myid==0) printf("before smoothit\n");
		smoothit(den,nx,ny,nz,esmooth);

		if(myid==0) printf("before stat3d\n");
		stat3d(den,nx,ny,nz,ismooth,0);

		if(myid==0) printf("before gradx\n");
		gradx(den,nx,ny,nz);

		if(myid==0) printf("before stat3d\n");
		stat3d(den,nx,ny,nz,ismooth,1);


		/*
		MPI_Finalize();return 0;
		FILEDUMP(std,stdg,outnum,ismooth);
		*/
		free(den);

		if(myid==0) printf("before copyden\n");
		den = copyden(iden, nx, ny,nz);

		if(myid==0) printf("before smoothid\n");
		smoothit(den,nx,ny,nz,esmooth);
		if(myid==0) printf("before statistics\n");
		Statistics(den, nx, ny, nz);


		if(myid==0) printf("before sampling\n");
		float *densample = Sampling(den,nx,ny,nz,sx,sy,sz);
		if(myid==0) {
			qsort(densample,npsample,sizeof(float),sorttype);
		}
		MPI_Bcast(densample,npsample,MPI_FLOAT,0,MPI_COMM_WORLD);

		if(myid==0){
			for(i=0;i<81;i++){
				char outsub[190]; sprintf(outsub,"subrodgenus%.2d.out.imock%.2d",ismooth,i);
				FILE *wwp = fopen(outsub,"w");
				fclose(wwp);
			}
		}


		int readflag = 1;
		do {
			if(myid==0) {
				int ncount;
				ncount = fscanf(fp,"%g %g %g\n",&label,&vmfrac,&delta);
				if(ncount !=3) readflag = 0;
				if(readflag) printf("%g %g %g\n",label,vmfrac,delta);
			}
			MPI_Bcast(&readflag,1,MPI_INT,0,MPI_COMM_WORLD);
			if(readflag){
				MPI_Bcast(&vmfrac,1,MPI_FLOAT,0,MPI_COMM_WORLD);
				long long iisort = (1.L-vmfrac)*npsample;
				if(iisort >= npsample) iisort = npsample-1;
				else if(iisort < 0) iisort = 0;
				rhoc = densample[iisort];
				if(myid==0) printf("rhoc = %g\n",rhoc);
				int *hilo = (int *)malloc(sizeof(int)*(nx)*(ny)*(local_nz+2));
				CalHighLow(den,nx,ny,nz,rhoc,hilo);
				if(pein == 'p'){
					paddinghighlow(hilo,nx,ny,nz);
				}
				MPI_Barrier(MPI_COMM_WORLD);

				summ(hilo,nx,ny,nz,ismooth);
				free(hilo);
			}
		} while(readflag);

		free(densample);free(den);
		if(myid==0){
			fclose(fp);
			fclose(wp1);
			fclose(wp2);
		}

	}
	}
	MPI_Finalize();
	return 0;
}

void summ(int *hilo, int nx, int ny, int nz, int ismooth){
	int d[2][256] = { 
		{ 
			0, 1, 1, 0, 1, 2, 0,-1, 1, 0,-2,-1, 0,-1,-1, 0, 1, 0, 2,-1, 
			2, 1, 1,-2, 2,-1,-1,-2, 1,-2, 0,-1, 1,-2, 0,-1, 2,-1,-1,-2, 
			2,-3,-3,-2, 1,-4,-2,-1, 0,-1,-1, 0, 1, 0,-2,-1, 1,-2,-4,-1, 
			0,-3,-3, 0, 1, 2,-2,-3, 0, 1,-1,-2,-2,-3,-1, 0,-1,-2,-2,-1, 
			2, 1,-1,-4, 1, 0, 0,-3,-1,-4, 0,-1, 0,-3,-1,-2, 0,-3,-1,-2, 
			-1,-4, 0,-1,-3,-8, 0, 1,-2,-7,-1, 0,-1,-2,-2,-1,-2,-3,-1, 0, 
			-4,-7,-1, 2,-3,-6,-2, 1, 1,-2,-6,-3, 2,-1,-7,-4, 0,-1,-3,-2, 
			-1,-2,-2,-1, 0,-1,-7,-2, 1, 0,-8,-3,-1, 0,-4,-1,-2,-1,-3, 0, 
			-2,-1,-3, 0,-1, 0,-4,-1,-3, 0, 0, 1,-4,-1, 1, 2,-1,-2,-2,-1, 
			0,-1,-3,-2,-2,-1, 1, 0,-3,-2, 2, 1, 0,-3,-3, 0,-1,-4,-2, 1, 
			-1,-2, 0, 1, 0,-1,-1, 0,-1,-2,-4, 1,-2,-3,-3, 2,-2,-1,-1, 2, 
			-1, 0,-2, 1,-1, 0,-2, 1,-2,-1,-1, 2,-2, 1, 1, 2,-1, 2, 0, 1, 
			0,-1,-1, 0,-1,-2, 0, 1,-1, 0, 2, 1, 0, 1, 1, 0
		},
		{ 
			0, 1, 1, 0, 1,-2, 0,-1, 1, 0, 2,-1, 0,-1,-1, 0, 1, 0, 2,-1, 
			2,-3, 1,-2, 2,-1, 3,-2, 1,-2, 0,-1, 1,-2, 0,-1, 2,-5,-1,-2, 
			2,-3, 1,-2, 1,-4,-2,-1, 0,-1,-1, 0, 1,-4,-2,-1, 1,-2, 0,-1, 
			0,-3,-3, 0, 1,-6,-2,-3, 0,-3,-1,-2,-2,-3,-5, 0,-1,-2,-2,-1, 
			-2,-3,-5, 0,-3, 0,-4, 1,-5, 0,-8, 3,-4, 1,-5, 2, 0,-3,-1,-2, 
			-1, 0, 0,-1,-3, 0,-4, 1,-2, 1,-1, 0,-1,-2,-2,-1,-2, 1,-1, 0, 
			-4, 1,-5, 2,-3, 2,-2, 1, 1,-2, 2,-3, 2,-5, 1,-4, 0,-1, 1,-2, 
			-1,-2,-2,-1, 0,-1, 1,-2, 1,-4, 0,-3,-1, 0, 0,-1,-2,-1,-3, 0, 
			2,-5, 1,-4, 3,-8, 0,-5, 1,-4, 0,-3, 0,-5,-3,-2,-1,-2,-2,-1, 
			0,-5,-3,-2,-2,-1,-3, 0,-3,-2,-6, 1, 0,-3,-3, 0,-1, 0,-2, 1, 
			-1,-2,-4, 1, 0,-1,-1, 0,-1,-2,-4, 1,-2, 1,-3, 2,-2,-1,-5, 2, 
			-1, 0,-2, 1,-1, 0,-2, 1,-2, 3,-1, 2,-2, 1,-3, 2,-1, 2, 0, 1, 
			0,-1,-1, 0,-1, 2, 0, 1,-1, 0,-2, 1, 0, 1, 1, 0
		}
	};
	int add = 0;
	int zadd = 0;
	if(myid == nid-1 && pein == 'i') zadd = 1;
	if(pein == 'i') add = 1;
	int i,j,k;
	int kk;
	long long *subdeficit = (long long*)calloc(nsubx*nsuby*nsubz,sizeof(long long));
	long long *tsubdeficit = (long long*)calloc(nsubx*nsuby*nsubz,sizeof(long long));
	float *subgenus = (float*)calloc(nsubx*nsuby*nsubz,sizeof(float));
	long long nvert,deficit,Lnvert, Ldeficit;
	Lnvert = Ldeficit = 0;
	for(k=0;k<local_nz-zadd;k++){
		for(j=0;j<ny-add;j++){
			int jp1 = (j+1+ny)%ny;
			int itype = (j+1)%2 ;
			for(i=0;i<nx-add;i++){
				int ip1 = (i+1+nx)%nx;
				int s1,s2;
				s1 = 8*hilo[i+nx*(j+ny*k)] + 4*hilo[i+nx*(j+ny*(k+1))] +
					2*hilo[i+nx*(jp1+ny*(k+1))] + hilo[i+nx*(jp1+ny*(k))];
				s2 = 8*hilo[ip1+nx*(j+ny*k)] + 4*hilo[ip1+nx*(j+ny*(k+1))] +
					2*hilo[ip1+nx*(jp1+ny*(k+1))] + hilo[ip1+nx*(jp1+ny*(k))];
				int iodd = (i+1)%2;
				int ieven = 1 - iodd;
				int s = s1*(iodd*16+ieven) + s2*(ieven*16+iodd);
				if( s>=0) {
					Lnvert ++;
					Ldeficit += d[itype][s];
					int ii = i/cubesizex;
					int jj = j/cubesizey;
					int kk = (k+local_z_start)/cubesizez;
					if(ii<nsubx && jj<nsuby &&kk<nsubz) subdeficit[ii+nsubx*(jj+nsuby*kk)] += d[itype][s];
				}
			}
		}
	}
	MPI_Reduce(&Ldeficit,&deficit,1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&Lnvert,&nvert,1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(&deficit,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(&nvert,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
	genus = -1.*deficit/8.;
	MPI_Reduce(subdeficit,tsubdeficit,nsubx*nsuby*nsubz,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(tsubdeficit,nsubx*nsuby*nsubz,MPI_LONG_LONG,0,MPI_COMM_WORLD);
	for(i=0;i<nsubx*nsuby*nsubz;i++) subgenus[i] = -1.*tsubdeficit[i]/8.;



	double vfrac = vol/vtot;
	double mfrac = mass/mtot;

	if(myid==0){
		fprintf(wp1,"%8.2f %14.4f %g %g %g\n",label,genus,rhoc,vfrac,mfrac);
		fprintf(wp2,"number of interior vertices = %ld\n",nvert);
		fprintf(wp2,"criterion density = %g\n",rhoc);
		fprintf(wp2,"high density fractional volume = %g\n",vfrac);
		fprintf(wp2,"high density fractional mass = %g\n",mfrac);
		fprintf(wp2,"genus of boundary surface = %g\n",genus);
		fprintf(stdout,"number of interior vertices = %ld\n",nvert);
		fprintf(stdout,"criterion density = %g\n",rhoc);
		fprintf(stdout,"high density fractional volume = %g\n",vfrac);
		fprintf(stdout,"high density fractional mass = %g\n",mfrac);
		fprintf(stdout,"genus of boundary surface = %g\n",genus);
		fflush(wp1); fflush(wp2);
		/*
		fwrite(&label,sizeof(float),1,wwp);
		fwrite(&rhoc,sizeof(float),1,wwp);
		fwrite(&vfrac,sizeof(float),1,wwp);
		fwrite(&mfrac,sizeof(float),1,wwp);
		int msub = nsub*nsub*nsub;
		fwrite(&msub,sizeof(int),1,wwp);
		fwrite(subgenus,sizeof(float),msub,wwp);
		*/
		for(i=0;i<81;i++){
			char outsub[190]; sprintf(outsub,"subrodgenus%.2d.out.imock%.2d",ismooth,i);
			FILE *wwp = fopen(outsub,"a");
			fprintf(wwp,"%8.2f %14.4f %g %g %g\n",label,subgenus[i],rhoc,vfrac,mfrac);
			fclose(wwp);
		}
	}
	free(subdeficit);free(tsubdeficit);free(subgenus);
}

void paddinghighlow(int *hilo,int nx,int ny, int nz){
	int i,j,k;
	long long sendsize, recvsize;
	MPI_Status status;
	sendsize = nx*ny;
	recvsize = nx*ny;
	int dest,src;
	int stag, rtag;
	stag = rtag = 0;
	int *arr_local_nz = (int * ) malloc(sizeof(int)*nid);
	if(myid==0){
		arr_local_nz[0]= local_nz;
		for(i=1;i<nid;i++) MPI_Recv(arr_local_nz+i,1,MPI_INT,i,rtag,MPI_COMM_WORLD,&status);
	}
	else MPI_Send(&local_nz,1,MPI_INT,0,stag,MPI_COMM_WORLD);
	MPI_Bcast(arr_local_nz,nid,MPI_INT,0,MPI_COMM_WORLD);
	/* sending data to upper rank */
	for(i=1;i<nid;i++){
		dest = (myid+i+nid)%nid;
		if(arr_local_nz[dest]>0){
			break;
		}
	}
	for(i=1;i<nid;i++){
		src = (myid-i+nid)%nid;
		if(arr_local_nz[src]>0){
			break;
		}
	}
	if(local_nz>0){
		MPI_Sendrecv(hilo+sendsize*(local_nz-1),sendsize,MPI_INT, dest, stag,
				hilo+recvsize*(local_nz+1),recvsize,MPI_INT, src, rtag,MPI_COMM_WORLD,&status);
	}
	/* sending data to lower rank */
	for(i=1;i<nid;i++){
		dest = (myid-i+nid)%nid;
		if(arr_local_nz[dest]>0){
			break;
		}
	}
	for(i=1;i<nid;i++){
		src = (myid+i+nid)%nid;
		if(arr_local_nz[src]>0){
			break;
		}
	}
	if(local_nz>0){
		MPI_Sendrecv(hilo,sendsize,MPI_INT, dest, stag,
				hilo+recvsize*(local_nz),recvsize,MPI_INT, src, rtag,MPI_COMM_WORLD,&status);
	}
	free(arr_local_nz);
}


void CalHighLow(float *den, int nx, int ny, int nz, float rhoc, int *hilo){
	int i,j,k;
	int nnx = 2*(nx/2+1);
	long long high = 0, thigh;
	double tmass;
	vol = mass = 0;
	for(k=0;k<local_nz;k++){
		for(j=0;j<ny;j++){
			for(i=0;i<nx;i++){
				if(den[i+nnx*(j+ny*k)]>rhoc){
					hilo[i+nx*(j+ny*k)] = 1;
					high ++;
					mass += den[i+nnx*(j+ny*k)];
				}
				else if(den[i+nnx*(j+ny*k)]>=0){
					hilo[i+nx*(j+ny*k)] = 0;
				}
				else {
					hilo[i+nx*(j+ny*k)] = -256;
				}
			}
		}
	}
	MPI_Reduce(&high,&thigh,1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&mass,&tmass,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	high = thigh;
	mass = tmass;
	MPI_Bcast(&high,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(&mass,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
	vol = high;
}


void shuffleitandsample(float *den, int nx, int ny, int nz,
		float *sample, int sx, int sy, int sz, float s_local_nz){
	int nnx = 2*(nx/2+1);
	float *buff = (float*)malloc(sizeof(float)*nx*ny*local_nz);
	int i,j,k;
	for(k=0;k<local_nz;k++){
		for(j=0;j<ny;j++){
			for(i=0;i<nx;i++){
				buff[i+nx*(j+ny*k)] = den[i+nnx*(j+ny*k)];
			}
		}
	}
	long long ii,np;
	np = nx*ny*local_nz;
	for(ii=0;ii<np;ii++){
		long long jj = np*drand48();
		float tmp;
		if(jj>=np) jj = np -1;
		if(jj<0) jj=0;
		tmp = buff[ii];
		buff[ii] = buff[jj];
		buff[jj] = tmp;
	}

	for(ii=0;ii<sx*sy*s_local_nz;ii++) sample[ii] = buff[ii];
	free(buff);
}

float *Sampling(float *den, int nx,int ny,int nz,int sx,int sy,int sz){
	int i;
	float *densample,s_local_nz;
	s_local_nz = local_nz*(double) sz/(double)nz;
//	densample = (float *)malloc(sizeof(float)*sx*sy*s_local_nz);
	densample = (float *)malloc(sizeof(float)*sx*sy*sz);
	shuffleitandsample(den,nx,ny,nz,densample,sx,sy,sz,s_local_nz);
	if(myid==0){
		/*
		densample = (float*)realloc(densample,sizeof(float)*sx*sy*sz);
		*/
		long long mp;
		npsample = sx*sy*s_local_nz;
		MPI_Status status;
		for(i=1;i<nid;i++){
			MPI_Recv(&mp,1,MPI_LONG_LONG,i,i,MPI_COMM_WORLD,&status);
			MPI_Recv(densample+npsample,mp,MPI_FLOAT,i,i,MPI_COMM_WORLD,&status);
			npsample += mp;
		}
	}
	else {
		long long mp = sx*sy*s_local_nz;
		MPI_Send(&mp,1,MPI_LONG_LONG,0,myid,MPI_COMM_WORLD);
		MPI_Send(densample,mp,MPI_FLOAT,0,myid,MPI_COMM_WORLD);
	}
	MPI_Bcast(&npsample,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
	if(myid!=0) densample = (float*)malloc(sizeof(float)*npsample);
	if(myid==0) densample = (float*)realloc(densample,sizeof(float)*npsample);
	return densample;
}


void FILEDUMP(double std, double stdg, int outnum, int ismooth){
	char outfile[190];
	FILE *wp;
	sprintf(outfile,"statistics%.2d.dat%.2d",outnum,ismooth);
	wp = fopen(outfile,"w");
	fprintf(wp,"%g %g\n",std,stdg);
	fclose(wp);
}


float *copyden(float *iden,int nx, int ny, int nz){
	int i,j,k;
	float *res;
	int nnx = 2*(nx/2+1);
	res = (float *)malloc(sizeof(float)*nnx*ny*local_nz);

	for(k=0;k<local_nz;k++){
		for(j=0;j<ny;j++){
			for(i=0;i<nx;i++){
				res[i+nnx*(j+ny*k)] = iden[i+nx*(j+ny*k)];
			}
		}
	}
	return res;
}

void smoothit(float *den,int nx,int ny,int nz,float esmooth){
	float lambda = esmooth/(boxsize/nx);
	lambda = lambda/nx;
	float expfac = -(M_PI*lambda)*(M_PI*lambda);
	int i,j,k;
	int nnx = 2*(nx/2+1);

	rfftwnd_mpi(plan,1,den,NULL,FFTW_TRANSPOSED_ORDER);

	for(j=0;j<local_ny_after_transpose;j++){
		float j1 = j + local_y_start_after_transpose;
		if(j1 > ny/2) j1 = ny-j1;
		for(k=0;k<nz;k++){
			float k1 = k ;
			if(k1>nz/2) k1 = nz-k1;
			for(i=0;i<nx/2+1;i++){
				float wave = i*i+k1*k1+j1*j1;
				float fac = exp(wave*expfac);
				int i2 = i *2;
				den[i2  +nnx*(k+nz*j)] *= fac;
				den[i2+1+nnx*(k+nz*j)] *= fac;
			}
		}

	}
	rfftwnd_mpi(iplan,1,den,NULL,FFTW_TRANSPOSED_ORDER);
	double rincube = 1.L/(double)nx/(double)ny/(double)nz;
	for(i=0;i<nnx*ny*local_nz;i++) den[i] *= rincube;
}


void stat3d(float *den, int nx, int ny, int nz,int ismooth,int iflag){
	int i,j,k;
	int nnx=2*(nx/2+1);
	int ii,jj,kk;
	double *submean,*substd,*tsubmean,*tsubstd;
	double a,s;
	long msub = nsubx*nsuby*nsubz;
	submean = (double*)calloc(msub,sizeof(double));
	substd = (double*)calloc(msub,sizeof(double));
	tsubmean = (double*)calloc(msub,sizeof(double));
	tsubstd = (double*)calloc(msub,sizeof(double));
	if(myid==0) printf(" %d %d %d %d : %d: %g %g %g\n",nx,nsubx,nsuby,nsubz,cubesizex,den[0],den[3456],den[23485]);
	for(k=0;k<local_nz;k++){
		kk = (k+local_z_start)/cubesizez;
		for(j=0;j<ny;j++) {
			jj = j/cubesizey;
			for(i=0;i<nx;i++){
				ii = i/cubesizex;
				if(ii<nsubx && jj<nsuby && kk<nsubz) {
					float d = den[i+nnx*(j+ny*k)];
					submean[ii+nsubx*(jj+nsuby*kk)] += d;
					substd[ii+nsubx*(jj+nsuby*kk)] += d*d;
				}

			}
		}
	}

	MPI_Reduce(submean,tsubmean,msub,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(substd,tsubstd,msub,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	if(myid==0){
		for(i=0;i<msub;i++){
			tsubmean[i] = tsubmean[i]/cubesizex/cubesizey/cubesizez;
			tsubstd[i] = tsubstd[i]/cubesizex/cubesizey/cubesizez - tsubmean[i]*tsubmean[i];
			tsubstd[i] = sqrt(tsubstd[i]);
		}
		for(i=0;i<81;i++){
			char outsub[190];
			if(iflag ==0) sprintf(outsub,"sigma0_subplategenus%.2d.out.imock%.2d",ismooth,i);
			else if(iflag ==1) sprintf(outsub,"sigma1_subplategenus%.2d.out.imock%.2d",ismooth,i);
			FILE *wwp = fopen(outsub,"w");
			fprintf(wwp,"%g %g\n",tsubmean[i],tsubstd[i]);
			fclose(wwp);
		}
	}
	free(submean);free(substd);free(tsubmean);free(tsubstd);
}


void gradx(float *den,int nx, int ny, int nz){
	int i,j,k;
	int nnx = 2*(nx/2+1);

	rfftwnd_mpi(plan,1,den,NULL,FFTW_TRANSPOSED_ORDER);
	for(j=0;j<local_ny_after_transpose;j++){
		float j1 = j + local_y_start_after_transpose;
		if(j1 > ny/2) j1 = j1-ny;
		for(k=0;k<nz;k++){
			float k1 = k ;
			if(k1>nz/2) k1 = k1 - nx;
			for(i=0;i<nx/2+1;i++){
				float wave = i*i+j1*j1+k1*k1;
				if(wave ==0) {
					den[2*i  +nnx*(k+nz*j)] = 0;
					den[2*i+1+nnx*(k+nz*j)] = 0;
				}
				else {
					float dreal = den[2*i+1+nnx*(k+nz*j)];
					den[2*i  +nnx*(k+nz*j)] = -i *den[2*i+1+nnx*(k+nz*j)];
					den[2*i+1+nnx*(k+nz*j)] = i*dreal;
				}
			}
		}
	}
	rfftwnd_mpi(iplan,1,den,NULL,FFTW_TRANSPOSED_ORDER);
	double pixel = boxsize/(double)nz;
	double rincube = 2.*M_PI/(double)nx/(double)ny/(double)nz/pixel/(double)nz;
	for(i=0;i<nnx*ny*local_nz;i++) den[i] *= rincube;
}


void Statistics(float *den, int nx, int ny, int nz){
	double mt,vt;
	mt = vt=0;
	int i,j,k;
	int nnx = 2*(nx/2+1);
	for(k=0;k<local_nz;k++){
		for(j=0;j<ny;j++){
			for(i=0;i<nx;i++){
				if(den[i+nnx*(j+ny*k)]>0) {
					vt += 1.;
					mt += den[i+nnx*(j+ny*k)];
				}
			}
		}
	}
	MPI_Reduce(&mt,&mtot,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&vt,&vtot,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(&mtot,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&vtot,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
}
float *meshdenread(char *filename, int *nx, int *ny, int *nz){
	float *den;
    int myid,nid;
    char inhead[190],infile[190];
    FILE *fp;
    int mx,my,mz;
    int itag = 1;
    MPI_Status status;
    long long nsize;
    int i,j,k;

    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&nid);
    if(myid == nid-1){
        printf("P%d is opening file %s\n",myid,filename);
        fp = fopen(filename,"r");
        fread(nx,sizeof(int),1,fp);
        fread(ny,sizeof(int),1,fp);
        fread(nz,sizeof(int),1,fp);
		printf("nx/ny/nz = %d %d %d\n",*nx,*ny,*nz);
    }
    MPI_Bcast(nx,1,MPI_INT,nid-1,MPI_COMM_WORLD);
    MPI_Bcast(ny,1,MPI_INT,nid-1,MPI_COMM_WORLD);
    MPI_Bcast(nz,1,MPI_INT,nid-1,MPI_COMM_WORLD);
	mx = *nx;
	my = *ny;
	mz = *nz;
	den = (float*)malloc(sizeof(float)*mx*my*local_nz);
    if(myid == nid-1){
        for(i=0;i<nid;i++){
            if(i != myid){
				float *localden;
				int localsize;
				MPI_Recv(&localsize,1,MPI_INT,i,0,MPI_COMM_WORLD,&status);
        		nsize = (long long)(mx)*(long long)(my)*(long long)localsize;
				localden = (float*)malloc(sizeof(float)*nsize);
				printf("P%d is reading %d slices for p%d\n",myid,localsize,i);
            	fread(localden,sizeof(float),nsize,fp);
                MPI_Send(localden,nsize,MPI_FLOAT,i,0,MPI_COMM_WORLD);
				free(localden);
            }
			else{
				nsize = (long long)mx*(long long)my*(long long)local_nz;
				fread(den,sizeof(float),nsize,fp);

			}
        }
        fclose(fp);

    }
    else {
        nsize = (long long)mx*(long long)my*(long long)local_nz;
        MPI_Send(&local_nz,1,MPI_INT,nid-1,0,MPI_COMM_WORLD);
        MPI_Recv(den,nsize,MPI_FLOAT,nid-1,0,MPI_COMM_WORLD,&status);
    }
    nsize = (long long)mx*(long long)my*(long long)local_nz;
    if(nsize>0){
        double mean,std;
        mean = std = 0;
        for(j=0;j<nsize;j++){
			/*
            den[j] += 5.;
			*/
            mean += den[j];
            std += den[j]*den[j];
        }
        mean = mean / nsize;
        std = std/nsize - mean*mean;
        std = sqrt(std);
        printf("P%d has mean/std = %g %g\n",myid,mean,std);
    }
	return den;
}
