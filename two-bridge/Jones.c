#include<stdio.h>
#include<complex.h>
#include<math.h>

int main(void) {
	int n;
	double _Complex q;
	printf("n�̒l����͂��Ă��������B");
	scanf("%d",&n);
	q = cexp(2*M_PI*I/n);
	
	printf("%.24lf\n",cabs(cpow(q,n/2)-cpow(q,(-1)*n/2)));
	
	return 0;
}