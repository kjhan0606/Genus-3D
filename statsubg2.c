#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>



int main(int argc, char **argv){
	long long i,j,k;
	FILE *fp;

	float *genus;
	int nnu=501;



	int ismooth=0; 
	int mp;
	char outfile[190];
	sprintf(outfile,"statsubgenus");
	FILE *wp = fopen(outfile,"w");

	{
		char infile[190];
		sprintf(infile,"subgenus%.2d.bin",ismooth);
		fp = fopen(infile,"r");
		float a;
		float label,rhoc,vfrac,mfrac;
		int nsub;
		while((mp=fread(&label,sizeof(float),1,fp)) >0){
			fread(&rhoc,sizeof(float),1,fp);
			fread(&vfrac,sizeof(float),1,fp);
			fread(&mfrac,sizeof(float),1,fp);
			fread(&nsub,sizeof(int),1,fp);
			genus = (float*)malloc(sizeof(float)*nsub);
			fread(genus,sizeof(float),nsub,fp);
			double mean,std;
			mean = std = 0;
			int isub,ncount=0;
			for(isub=0;isub<nsub;isub+=27){
				mean += genus[isub];
				std += genus[isub] * genus[isub];
				ncount ++;
			}
			mean = mean/ncount;
			std = std/ncount - mean*mean;
			std = sqrt(std);
			fprintf(wp,"%10.4f %14.4f %14.4f\n",label,mean,std);

			free(genus);
		}

		fclose(fp);
	}
	fclose(wp);
}
