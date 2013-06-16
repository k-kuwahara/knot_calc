#include<stdio.h>
#include<complex.h>
#include<math.h>

/*ーーーープロトタイプ宣言ーーーー*/
double _Complex Q_minus(int N, int n, double _Complex q);
double _Complex Q_plus(int N, int n, double _Complex q);


/*ーーーーメイン関数ーーーー*/
int main(void){
	int N,n;
	for(N=2162; N<2164; N++){
		double _Complex q = cexp(2*M_PI*I/N);
		double _Complex a = cexp(M_PI*I/N);
		long double _Complex J=0;
		
		for(n=0; n<=N-1; n++){
			J = J + pow(-1,n) * cpow(a,-n*(n+1)) * Q_minus(N,n,q) * Q_plus(N,n,q);
		}
		printf("N=%d：%.60lf\n\n",N,log(cabs(J))/N);
	}
	return 0;
}

/*ーーーーサブ関数ーーーー*/
double _Complex Q_minus(int N, int n, double _Complex q) {
	int i,l;			//カウント用
	long double _Complex re_value = 1;	//戻り値用の変数宣言
	if(n == 0) return 1;
	else {
		for(l=1; l<=n; l++){
			re_value *= (1 - cpow(q,1-N)*cpow(q,n-l));
		}
	}
	printf("Q_minus:%lf  %lf\n",creal(re_value),cimag(re_value));
	return re_value;
}

double _Complex Q_plus(int N, int n, double _Complex q) {
	int i,l;			//カウント用
	long double _Complex re_value = 1;	//戻り値用の変数宣言
	if(n == 0) return 1;
	else {
		for(l=1; l<=n; l++){
			re_value *= (1 - cpow(q,1+N)*cpow(q,n-l));
		}
	}
	printf("Q_plus:%lf  %lf\n\n",creal(re_value),cimag(re_value));
	return re_value;
}