#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

struct input {int Ntot; double Nin;};

double randNumber(void) {unsigned int seed; double output = 1.0*rand_r(&seed)/RAND_MAX; return output;}

void* bar(void* arg) {
	struct input *p = (struct input*)arg;
	int N = (*p).Ntot;
	int count = 0;
	for(int i=0; i<N; i++) {
		double x = randNumber(); double y = randNumber();
		if (pow(x,2)+pow(y,2)<1) {
			count++;
		}
	}
	p->Nin=count;
return NULL;
}

int main() {
	int Ntot = 1e8;
	double sum1 = 0, sum2 = 0, sum3 = 0;
	pthread_t t1,t2,t3;
	pthread_attr_t* attributes = NULL;
	struct input p1 = {.Ntot = Ntot,.Nin = sum1};
	struct input p2 = {.Ntot = Ntot,.Nin = sum2};
	struct input p3 = {.Ntot = Ntot,.Nin = sum3};
	pthread_create(&t1, attributes, bar, (void*)& p1);
	pthread_create(&t2, attributes, bar, (void*)& p2);
	pthread_create(&t3, attributes, bar, (void*)& p3);
	void* returnvalue = NULL;
	pthread_join(t1,returnvalue);
	pthread_join(t2,returnvalue);
	pthread_join(t3,returnvalue);
	double NinAVG = (p1.Nin + p2.Nin + p3.Nin) / 3;
	printf("A: Pi value is %10g\n",4*NinAVG/Ntot);
return 0;
}
