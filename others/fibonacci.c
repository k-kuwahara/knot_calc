#include<stdio.h>

int main(void){
	int i, n;
	long double  prev=0, next=1, buf=0;

	printf("type one natural number\t");
	scanf("%d", &n);

	for(i=2; i<n; i++){
		buf = next + prev;
		prev = next;
		next = buf;

		printf("%d: %.0Lf\n", i, buf);
	}

	return 0;
    
    
//23601まで計算可能
}