#include<stdio.h>

int main(void){
	int i,n,sosu=2001,answer=1;
	
//	printf("自然数を一つ入力して下さい。");
	printf("type one natural number\n");
	scanf("%d",&n);
	printf("\n");
	
	if(n <= 2){
		printf("mistype!!\n");
	}
	else if(n<sosu){
		printf("The number is lower than a prime number!!\n");
	}
	if(n > 2){
		do{
			/*奇素数の判定*/
			for(i = 3; i * i <= sosu; i += 2){
//				printf("%d:%d  ",i, sosu);
				if(sosu % i == 0){
//					printf("I'm here!!\n");
					answer=0;
					break;
				}
			}
//			printf("\n%d\n",answer);
			if(answer == 1)	printf("%d ",sosu);
			sosu += 2;
			answer = 1;
		}while(sosu<=n);
	}
	printf("\n");
	return 0;
}