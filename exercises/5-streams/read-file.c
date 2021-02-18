#include<stdio.h>
#include<math.h>

int main(int argc, char **argv){
	if(argc<2) fprintf(stderr,"%s: there were no arguments\n",argv[0]);
	if(argc>3) fprintf(stderr,"%s: there were too many arguments",argv[0]);
	else {
		FILE* inputfile = fopen(argv[1],"r");
		FILE* outputfile = fopen(argv[2],"w");
		double x;
		int items;
		items = fscanf(inputfile, "%lg",&x);
		while (items!=EOF) {
			fprintf(outputfile, "x = %lg, sin(x) = %lg, cos(x) = %lg\n", x, sin(x), cos(x));
			items = fscanf(inputfile, "%lg",&x);
		}
	}
return 0;
}
