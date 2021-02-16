#include<math.h>

int equals(double a, double b, double tau, double epsilon){
	if (fabs(a-b)<tau) return 1;
	if (fabs(a-b)/(fabs(a)+fabs(b))) return 1;
	else return 0;}
