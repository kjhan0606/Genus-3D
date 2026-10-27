void my_MPI_Send(void *, long , MPI_Datatype , int , int , MPI_Comm );
void my_MPI_Recv(void *, long , MPI_Datatype , int , int , MPI_Comm , MPI_Status *);
void my_MPI_Sendrecv(void *, long, MPI_Datatype , int , int , void *, long , 
		MPI_Datatype , int , int , MPI_Comm , MPI_Status *);

typedef struct PrimeNumber{
	int prime, factor;
} PrimeNumber;

#define min(a,b) ((a)<(b)? (a):(b))
#define max(a,b) ((a)>(b)? (a):(b))


#define ISUNSIGNED(a) (((typeof(a))-1) >= 0)


#ifdef __MAIN_PQSORT__

#define GenerateParallelComp(TYPE,subtype,mem)  \
	int Psort_##TYPE##_##subtype##_##mem(const void *a, const void *b){\
		TYPE *aa = (TYPE *)a;\
		TYPE *bb = (TYPE *)b;\
		if(aa->mem < bb->mem) return -1;\
		else if(aa->mem > bb->mem) return 1;\
		else return 0;\
	}\
	int Psort_##TYPE##_##subtype##_##mem##_L(const void *a, const void *lb, const void *lc){\
		if(((TYPE*)lb)->mem == ((TYPE*)lc)->mem) return 1;\
		else if(Psort_##TYPE##_##subtype##_##mem(a,lc) >=0 || Psort_##TYPE##_##subtype##_##mem(a,lb) ==-1) return 1;\
		else return 0;\
	}\
	int Psort_##TYPE##_##subtype##_##mem##_R(const void *a, const void *lb, const void *lc, int nid){\
		if(((TYPE*)lb)->mem == ((TYPE*)lc)->mem) return 1;\
		else if(Psort_##TYPE##_##subtype##_##mem(a,lc) >=0 || Psort_##TYPE##_##subtype##_##mem(a,lb) ==-1) {\
			TYPE *aa = (TYPE*) a;\
			TYPE *lbb = (TYPE*) lb;\
			TYPE *lcc = (TYPE*) lc;\
			int rdist;\
			rdist = (aa->mem - lbb->mem)/(lcc->mem-lbb->mem);\
			rdist = (rdist + nid)%nid;\
			if(rdist <= nid/2) return 1;\
			else return 0;\
		}\
		else return 0;\
	} 
#endif

void *onewaypqsort(void *, size_t *, size_t, int (*)(const void *, const void *, const void *, int ), 
		int (*)(const void *, const void *, const void *), 
		int (*)(const void *, const void *), void *, void *, MPI_Comm);
void *twowaypqsort(void *, size_t *, size_t, int (*)(const void *, const void *, const void *, int ), 
		int (*)(const void *, const void *, const void *), 
		int (*)(const void *, const void *), void *, void *, MPI_Comm);
void mpiqsort(void **, size_t *, size_t, int (*)(const void *, const void *), MPI_Comm);
void mpiqsort4noteven(void **,size_t *, size_t , int (*)(const void *, const void *), 
		float (*)(const void *, const void *), void (*)(const void *, const void *, float ), MPI_Comm );
