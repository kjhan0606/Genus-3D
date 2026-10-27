/* Version 0.5 of parallel sorting */
/* Copyright to Juhan Kim (kjhan@kias.re.kr) */
/* mpiqsort() */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<limits.h>
#include<mpi.h>
#include "mpiqsort.h"


#define SWAP(a,b,stmp) do{\
	memcpy(stmp,a,n_size);\
	memcpy(a,b,n_size);\
	memcpy(b,stmp,n_size);\
}while(0)

int getNextPrimeNumber(int number){
	int half = number/2;
	int i;
	/*
	for(i=half;i>=2;i--){
		if(number%i ==0) return i;
	}
	*/
	for(i=2;i<=half;i++){
		if(number%i ==0) return i;
	}
	return number;
}

int getprimenumber(int number,PrimeNumber *prime){
	int i,j,k,denom,nprime;
	int res;


	nprime = 0;
	prime[nprime].factor = 0;
	while((number %2)==0){
		prime[nprime].prime = 2;
		prime[nprime].factor ++;
		number = number / 2;
	}
	if(prime[nprime].factor >0) nprime ++;
	prime[nprime].factor = 0;
	for(i=3;i<=number;i+=2){
		while((number % i) ==0){
			prime[nprime].prime = i;
			prime[nprime].factor ++;
			number = number/i;
		}
		if(prime[nprime].factor > 0) nprime ++;
		prime[nprime].factor = 0;
	}
	return nprime;
}


void my_MPI_Sendrecv(void *sbase, long smem, MPI_Datatype datatype1, int dest, int stag,
	void *rbase, long rmem, MPI_Datatype datatype2, int src, int rtag, MPI_Comm comm, MPI_Status *status){
	long nmax, nchunk = INT_MAX/2, max;
	long i;
	MPI_Reduce(&smem, &nmax, 1,MPI_LONG, MPI_MAX, 0, comm);
	max = nmax;
	MPI_Reduce(&rmem, &nmax, 1,MPI_LONG, MPI_MAX, 0, comm);
	if(nmax > max) max = nmax;
	MPI_Bcast(&max, 1, MPI_LONG,0,comm);
	for(i=0;i<max;i+=nchunk){
		int ismem, irmem;
		long lsmem, lrmem;
		lsmem = smem -i;
		lrmem = rmem -i;
		if(lsmem > nchunk) ismem = nchunk;
		else if(lsmem <=0 ) ismem = 0;
		else ismem = lsmem;
		if(lrmem > nchunk) irmem = nchunk;
		else if(lrmem <=0 ) irmem = 0;
		else irmem = lrmem;
		MPI_Sendrecv((char*)sbase+i, ismem, datatype1,dest,stag, 
			(char*)rbase+i,irmem, datatype2, src, rtag, comm, status);
	}
}
void my_MPI_Send(void *sbase, long smem, MPI_Datatype datatype, int dest, int stag, MPI_Comm comm){
	long nchunk = INT_MAX/2;
	long i;
	for(i=0;i<smem;i+=nchunk){
		int ismem;
		long lsmem;
		lsmem = smem - i;
		if(lsmem > nchunk) ismem = nchunk;
		else if (lsmem <=0) ismem = 0;
		else ismem = lsmem;
		MPI_Send((char*) sbase+i, ismem, datatype, dest, stag, comm);
	}
}
void my_MPI_Recv(void *rbase, long smem, MPI_Datatype datatype, int src, int rtag, MPI_Comm comm, MPI_Status *status){
	long nchunk = INT_MAX/2;
	long i;
	for(i=0;i<smem;i+=nchunk){
		int ismem;
		long lsmem;
		lsmem = smem - i;
		if(lsmem > nchunk) ismem = nchunk;
		else if (lsmem <=0) ismem = 0;
		else ismem = lsmem;
		MPI_Recv((char*) rbase+i, ismem, datatype, src, rtag, comm, status);
	}
}

int GetRandom3(void *base, size_t nmem, size_t n_size, int (*compare)(const void *, const void *), void *pivot){
	int p1,p2,p3;
	char *a,*b,*c;
	if(nmem > 3) {
		p1 = (int) (nmem * (rand() / (RAND_MAX + 1.0)));
		p2 = (int) (nmem * (rand() / (RAND_MAX + 1.0)));
		while(p1 == p2) {p2 = (int) (nmem * (rand() / (RAND_MAX + 1.0)));}
		p3 = (int) (nmem * (rand() / (RAND_MAX + 1.0)));
		while(p3 == p1 || p3 == p2){ p3 = (int) (nmem * (rand() / (RAND_MAX + 1.0)));}
		a = (char*) base + p1*n_size;
		b = (char*) base + p2*n_size;
		c = (char*) base + p3*n_size;
		memmove((char*)pivot,a,n_size);
		memmove((char*)pivot+n_size,b,n_size);
		memmove((char*)pivot+2*n_size,c,n_size);
		return 3;
	}
	else if(nmem ==3){
		memmove((char*)pivot, (char*)base, 3*n_size);
		return 3;
	}
	else if(nmem ==2){
		memmove((char*)pivot, (char*)base, 2*n_size);
		return 2;
	}
	else if(nmem == 1){
		memmove((char*)pivot, (char*)base, n_size);
		return 1;
	}
	return 0;
}
int FindPivots(void *base, size_t nmem, size_t n_size, int (*comp)(const void *, const void *),
		int nsplit, void *pivot, int myid, int nid, MPI_Comm comm){
	int *displs = (int *)malloc(sizeof(int)*nid);
	int *ncounts = (int *)malloc(sizeof(int)*nid);
	void *temp = (void *)malloc(n_size*3);
	int npivot;

	int ipivot = GetRandom3(base, nmem, n_size, comp, temp);
	MPI_Gather(&ipivot,1, MPI_INT, ncounts, 1, MPI_INT, 0, comm);
	if(myid==0) {
		int i;
		for(i=0;i<nid;i++) ncounts[i] = ncounts[i] * n_size;
		displs[0] = 0;
		for(i=1;i<nid;i++) displs[i] = displs[i-1] + ncounts[i-1];
	}
	MPI_Gatherv(temp,ipivot*n_size, MPI_BYTE,pivot,ncounts, displs, MPI_BYTE,0,comm);
	if(myid==0) {
		int narray = (displs[nid-1]+ncounts[nid-1])/n_size;
		qsort(pivot, narray, n_size, comp);
		int istep = (narray)/(nsplit-1);
		int i,j;
		j = 0;
		for(i=istep/2;i<narray;i+=istep){
			memmove((char*)pivot + j*n_size, (char*)pivot+i*n_size,n_size);
			j++;
		}
		npivot = j;
	}
	MPI_Bcast(&npivot,1,MPI_INT,0,comm);
	MPI_Bcast(pivot,npivot*n_size,MPI_BYTE,0,comm);
	free(displs);
	free(ncounts);free(temp);
	return npivot;
}

void  determine_pivot(void *base, size_t nmem, size_t n_size, int (*compare)(const void *, const void *), void *pivot){
	char *a,*b,*c;
	int p1,p2,p3;
	if(nmem > 3) {
		p1 = (int) (nmem * (rand() / (RAND_MAX + 1.0)));
		p2 = (int) (nmem * (rand() / (RAND_MAX + 1.0)));
		while(p1 == p2) {p2 = (int) (nmem * (rand() / (RAND_MAX + 1.0)));}
		p3 = (int) (nmem * (rand() / (RAND_MAX + 1.0)));
		while(p3 == p1 || p3 == p2){ p3 = (int) (nmem * (rand() / (RAND_MAX + 1.0)));}
		a = (char*) base + p1*n_size;
		b = (char*) base + p2*n_size;
		c = (char*) base + p3*n_size;
	}
	else if(nmem ==3){
		a = (char *) base;
		b = a + n_size;
		c = b + n_size;
	}
	else if(nmem ==2){
		a = (char*)base;
		b = (char*)base + n_size;
		if(compare(a,b) >0) memmove(pivot, b, n_size);
		else memmove(pivot,a,n_size) ;
		return;
	}

	if(compare(a,b)>0){
		if(compare(b,c)>=0) memmove(pivot,b,n_size);
		else if(compare(a,c)>0) memmove(pivot,c,n_size);
		else if(compare(a,c)==0) memmove(pivot,b,n_size);
		else memmove(pivot,a,n_size);
	}
	else if(compare(a,b)==0){
		if(compare(a,c)>0) memmove(pivot,c,n_size);
		else memmove(pivot,a,n_size);
	}
	else{
		if(compare(a,c)>0) memmove(pivot,a,n_size);
		else if(compare(b,c)>0) memmove(pivot,c,n_size);
		else if(compare(b,c)==0) memmove(pivot,a,n_size);
		else memmove(pivot,b,n_size);
	}
}
void mpiqsort4noteven(void **ibase,size_t *mmem, size_t n_size, int (*Comp)(const void *, const void *), 
		float (*dist)(const void *, const void *), void (*SetVal)(const void *, const void *, float ), MPI_Comm Comm){
	MPI_Status status;
	int src,dest;
	long nrecv,nsend;
	long i,j,k;
	int myid,nid;
	size_t nmem = *mmem;
	void *base = *ibase;
	char swaptmp[n_size];
	int igroup,subgroupsize,subgroupid;
	size_t nowmem=0;
	void *rbase = NULL;
	char tlocalmin[n_size];
	char tlocalmax[n_size];

	MPI_Comm_rank(Comm,&myid);
	MPI_Comm_size(Comm,&nid);


	int nsubgroup = getNextPrimeNumber(nid);
	if(nsubgroup == 1) {
		qsort(*ibase,*mmem,n_size,Comp);
		return;
	}

	char tmin[n_size], tmax[n_size];
	char *left,*right;
	if(nmem >0) {
		memmove(tmin, (char*)base, n_size);
		memmove(tmax, (char*)base, n_size);
	}
	left = (char*)base;
	right = (char*)base + nmem*n_size;
	for(;left<right;left += n_size){
		if(Comp(left,tmax) >0) memmove(tmax,left,n_size);
		else if(Comp(left,tmin) <0) memmove(tmin,left,n_size);
	}
	char lmin[n_size],lmax[n_size];
	if(myid==0){
		for(i=1;i<nid;i++){
			MPI_Recv(lmin, n_size, MPI_BYTE, i,i, Comm, &status);
			MPI_Recv(lmax, n_size, MPI_BYTE, i,i, Comm, &status);
			if(Comp(lmin,tmin)<0) memmove(tmin,lmin, n_size);
			if(Comp(lmax,tmax)>0) memmove(tmax,lmax, n_size);
		}
	}
	else {
		MPI_Send(tmin, n_size, MPI_BYTE, 0,myid,Comm);
		MPI_Send(tmax, n_size, MPI_BYTE, 0,myid,Comm);
	}
	MPI_Bcast(tmin, n_size, MPI_BYTE, 0, Comm);
	MPI_Bcast(tmax, n_size, MPI_BYTE, 0, Comm);
	long lmem, tlmem;
	lmem = nmem;
	MPI_Reduce(&lmem,&tlmem,1, MPI_LONG, MPI_SUM,0,Comm);
	MPI_Bcast(&tlmem,1,MPI_LONG,0,Comm);
	long nbin = 16+nsubgroup+tlmem/4096/nid;
	long hist[nbin],thist[nbin];
	for(i=0;i<nbin;i++) hist[i] = 0;
	left = (char*)base;
	right = (char*)base + nmem*n_size;
	float stepsize = dist(tmax,tmin)/nbin;
#ifdef DEBUG
	if(myid==0) fprintf(stderr,"number of sampling bins %ld stepsize %g\n",nbin,stepsize);
#endif
	for(;left<right;left+=n_size){
		long ibin = dist(left, tmin)/stepsize;
		hist[ibin] ++;
	}
	MPI_Reduce(hist,thist, nbin,MPI_LONG,MPI_SUM,0,Comm);
	MPI_Bcast(thist,nbin,MPI_LONG,0,Comm);
	size_t tnmem=0;
	for(i=0;i<nbin;i++) tnmem += thist[i];
	size_t lnmem = (tnmem+nsubgroup-1)/nsubgroup;
	int id=0; long isum=0;
	char pivot[n_size*nid];
	for(i=0;i<nbin;i++){
		isum += thist[i];
		if(isum >= lnmem) {
			isum -= lnmem;
			SetVal(pivot+id*n_size, tmin, stepsize*i);
#ifdef DEBUG
			if(myid==0) fprintf(stderr,"pivot %d : %g nsubgroup= %d\n",
					id, *((float*)(pivot+id*n_size)),nsubgroup);
#endif
			id ++;
		}
	}

	subgroupsize = nid/nsubgroup;
	subgroupid = myid/subgroupsize;

	for(igroup=1;igroup<nsubgroup;igroup++){
		dest = (myid + subgroupsize*igroup + nid)%nid;
		src =  (myid - subgroupsize*igroup + nid)%nid;
		char localmin[n_size], localmax[n_size];
		int tgtsubgroupid = (subgroupid + igroup+nsubgroup)%nsubgroup;

		left = (char*)base;
		right = (char*)base + nmem*n_size;

		if(tgtsubgroupid==0){
			memmove(tlocalmax, (char*)pivot+tgtsubgroupid*n_size,n_size);
			for(;left<right;){
				if(Comp(left, tlocalmax)<0){
					right -= n_size;
					SWAP(left,right,swaptmp);
				}
				else left+=n_size;
			}
		}
		else if(tgtsubgroupid == nsubgroup-1){
			memmove(tlocalmin, (char*)pivot+(tgtsubgroupid-1)*n_size,n_size);
			for(;left<right;){
				if(Comp(left, tlocalmin) >= 0){
					right -= n_size;
					SWAP(left,right,swaptmp);
				}
				else left+=n_size;
			}
		}
		else {
			memmove(tlocalmin, (char*)pivot+(tgtsubgroupid-1)*n_size,n_size);
			memmove(tlocalmax, (char*)pivot+(tgtsubgroupid+0)*n_size,n_size);
			{
				for(;left<right;){
					if(Comp(left,tlocalmin)>=0 && Comp(left,tlocalmax)<0){
						right -= n_size;
						SWAP(left,right,swaptmp);
					}
					else left += n_size;
				}
			}
		}
		long  nsend = nmem-((char*)right-(char *)base)/n_size;
		MPI_Sendrecv(&nsend,1,MPI_LONG,dest,0,&nrecv,1,MPI_LONG,src,0,Comm,&status);
		if(rbase == NULL) rbase = (void *)malloc(sizeof(char)*nrecv*n_size);
		else rbase = (void *)realloc(rbase,sizeof(char)*(nrecv + nowmem) *n_size);
		my_MPI_Sendrecv(right,nsend*sizeof(char)*n_size,MPI_BYTE,dest,0,
				(char*)rbase+nowmem*n_size,nrecv*sizeof(char)*n_size,MPI_BYTE,src,0,Comm,&status);
		nmem = nmem-nsend;
		base = (void*)realloc(base, nmem * n_size);
		nowmem += nrecv;
		if(myid==0) fprintf(stderr,"+: %d : (%ld) <-> (%ld) %ld\n",igroup,nrecv, nsend,(long)nowmem); fflush(stdout);
	}
    rbase = (void *)realloc(rbase, sizeof(char)*(nowmem+nmem)*n_size);
    memmove((char*)rbase + nowmem*n_size, base, nmem*n_size);
    nowmem += nmem;
	free(base);
	{
		MPI_Comm newcom;
		int key = myid % subgroupsize;
		MPI_Comm_split(Comm, subgroupid, key, &newcom);
		mpiqsort4noteven(&rbase, &nowmem, n_size, Comp,dist,SetVal,newcom);
		MPI_Comm_free(&newcom);
	}
	*mmem = nowmem;
	*ibase =  rbase;
	return;
} 


void mpiqsort(void **ibase,size_t *mmem, size_t n_size, int (*Comp)(const void *, const void *), MPI_Comm Comm){
	MPI_Status status;
	int src,dest;
	long nrecv,nsend;
	long i,j,k;
	int myid,nid;
	size_t nmem = *mmem;
	void *base = *ibase;
	char swaptmp[n_size];
	int igroup,subgroupsize,subgroupid;
	size_t nowmem=0;
	void *rbase = NULL;
	char tlocalmin[n_size];
	char tlocalmax[n_size];

	MPI_Comm_rank(Comm,&myid);
	MPI_Comm_size(Comm,&nid);




	int nsubgroup = getNextPrimeNumber(nid);
	if(nsubgroup == 1) {
		qsort(*ibase,*mmem,n_size,Comp);
		return;
	}
	void *pivot = (void*)malloc(n_size*nid*3);
	int npivot = FindPivots(base, nmem, n_size, Comp,nsubgroup,pivot,myid, nid, Comm);
	subgroupsize = nid/nsubgroup;
	subgroupid = myid/subgroupsize;



	for(igroup=1;igroup<nsubgroup;igroup++){
		dest = (myid + subgroupsize*igroup + nid)%nid;
		src =  (myid - subgroupsize*igroup + nid)%nid;
		char localmin[n_size], localmax[n_size];
		int tgtsubgroupid = (subgroupid + igroup+nsubgroup)%nsubgroup;

		char *left,*right;
		left = (char*)base;
		right = (char*)base + nmem*n_size;

		if(tgtsubgroupid==0){
			memmove(tlocalmax, (char*)pivot+tgtsubgroupid*n_size,n_size);
			for(;left<right;){
				if(Comp(left, tlocalmax)<0){
					right -= n_size;
					SWAP(left,right,swaptmp);
				}
				else left+=n_size;
			}
		}
		else if(tgtsubgroupid == nsubgroup-1){
			memmove(tlocalmin, (char*)pivot+(tgtsubgroupid-1)*n_size,n_size);
			for(;left<right;){
				if(Comp(left, tlocalmin) >= 0){
					right -= n_size;
					SWAP(left,right,swaptmp);
				}
				else left+=n_size;
			}
		}
		else {
			memmove(tlocalmin, (char*)pivot+(tgtsubgroupid-1)*n_size,n_size);
			memmove(tlocalmax, (char*)pivot+(tgtsubgroupid+0)*n_size,n_size);
			{
				for(;left<right;){
					if(Comp(left,tlocalmin)>=0 && Comp(left,tlocalmax)<0){
						right -= n_size;
						SWAP(left,right,swaptmp);
					}
					else left += n_size;
				}
			}
		}
		long  nsend = nmem-((char*)right-(char *)base)/n_size;
		MPI_Sendrecv(&nsend,1,MPI_LONG,dest,0,&nrecv,1,MPI_LONG,src,0,Comm,&status);
		if(rbase == NULL) rbase = (void *)malloc(sizeof(char)*nrecv*n_size);
		else rbase = (void *)realloc(rbase,sizeof(char)*(nrecv + nowmem) *n_size);
		my_MPI_Sendrecv(right,nsend*sizeof(char)*n_size,MPI_BYTE,dest,0,
				(char*)rbase+nowmem*n_size,nrecv*sizeof(char)*n_size,MPI_BYTE,src,0,Comm,&status);
		nmem = nmem-nsend;
		base = (void*)realloc(base, nmem * n_size);
		nowmem += nrecv;
		if(myid==0) fprintf(stderr,"+: %d : (%ld) <-> (%ld) %ld\n",igroup,nrecv, nsend,(long)nowmem); fflush(stdout);
	}
    rbase = (void *)realloc(rbase, sizeof(char)*(nowmem+nmem)*n_size);
    memmove((char*)rbase + nowmem*n_size, base, nmem*n_size);
    nowmem += nmem;
	free(pivot);
	free(base);
	{
		MPI_Comm newcom;
		int key = myid % subgroupsize;
		MPI_Comm_split(Comm, subgroupid, key, &newcom);
		mpiqsort(&rbase, &nowmem, n_size, Comp,newcom);
		MPI_Comm_free(&newcom);
	}
	*mmem = nowmem;
	*ibase =  rbase;
	return;
} 

