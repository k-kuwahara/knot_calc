#include<stdio.h>
#include<complex.h>
#include<math.h>

/*ーーーグローバル変数ーーー*/
int m,N,P;
int k[4]={0,0,0,0};

/*ーーーサブ関数のプロトタイプ宣言ーーー*/
double _Complex f(int p, int n, double _Complex q);
double _Complex Q_minus(int n, double _Complex q);
double _Complex Q_plus(int n, double _Complex q);

double _Complex twist2(int len, double _Complex a, int n);
double _Complex twist3(int len, double _Complex a, int n);
double _Complex twist4(int len, double _Complex a, int n);

double _Complex comb(double _Complex a, int n);


/*ーーーメイン関数ーーー*/
int main(void) {
	int p,n;
	double _Complex q;
	long double _Complex Jones_sum=0;

	printf("Please type numbers of p(2~4)\n");
	scanf("%d",&p);
	P = p;
	
	for(N=2; N<1500; N++) {
		Jones_sum=1;
		q = cexp(2*M_PI*I/N);
		
		for(n=1; n<N; n++) {
			m=n;
			Jones_sum += f(p,n,q) * Q_minus(n,q) * Q_plus(n,q);
		}
//		printf("j_all:%.50lf\n",cabs(Jones_sum));
		printf("N=%d Jones:%.50lf\n",N,log(cabs(Jones_sum))/N);
	}
	return 0;
}


/*ーーーサブ関数ーーー*/
double _Complex f(int p, int n, double _Complex q) {
	double _Complex return_f=0, a = cexp(M_PI*I/N);
	
	if(p==2) {
		return_f = cpow(q,n) * twist2(p,a,n);
	}
	else if(p==3) {
		return_f = cpow(q,n) * twist3(p,a,n);
	}
	else if(p==4) {
		return_f = cpow(q,n) * twist4(p,a,n);
	}
	else {
		printf("Error:Number 'p' is wrong!\n");
		return_f = 0;
	}
	
	return return_f;
}

/*ーーーtwist2()ーーー*/
double _Complex twist2(int len, double _Complex a, int n) {
	int cnt,i;
	double _Complex re_value = 0;
	
	if(len==1) {
		k[len-1] = n;
//		for(i=0; i<2; i++) printf("%d ",k[i]);
//		printf(" end\n");
//		printf("n:%d\n",n);
		
		return cpow(a,(2*pow(k[1],2)+k[1]*k[0]+2*k[1]))*comb(a,m);
	}
	
	for(cnt=n;cnt>=0;cnt--){
		k[len-1] = cnt;
		re_value += twist2(len-1, a, n-cnt);
	}
	
	return re_value;
}

/*ーーーtwist3()ーーー*/
double _Complex twist3(int len, double _Complex a, int n) {
	int cnt,i;
	double _Complex re_value = 0;
	
	if(len==1) {
		k[len-1] = n;
//		for(i=0; i<3; i++) printf("%d ",k[i]);
//		printf(" end\n");
		
		return cpow(a,(2*pow(k[1],2)+4*pow(k[2],2)+k[1]*k[0]+5*k[1]*k[2]+k[2]*k[0]+2*k[1]+4*k[2]))*comb(a,m);
	}
	
	for(cnt=n;cnt>=0;cnt--){
		k[len-1] = cnt;
		re_value += twist3(len-1, a, n-cnt);
	}
	return re_value;
}

/*ーーーーtwist4()ーーーー*/
double _Complex twist4(int len, double _Complex a, int n) {
	int cnt,i;
	double _Complex re_value = 0;
	
	if(len==1) {
		k[len-1] = n;
//		for(i=0; i<4; i++) printf("%d ",k[i]);
//		printf(" end\n");
		
		return cpow(a,(6*pow(k[3],2)+4*pow(k[2],2)+2*pow(k[1],2)+9*k[3]*k[2]+5*k[3]*k[1]+k[3]*k[0]+
				5*k[2]*k[1]+k[2]*k[0]+k[1]*k[0]+2*k[1]+4*k[2]+6*k[3]))*comb(a,m);
	}
	
	for(cnt=n;cnt>=0;cnt--){
		k[len-1] = cnt;
		re_value += twist4(len-1, a, n-cnt);
	}
	return;	
}

/*ーーーーcombinationーーーー*/
double _Complex comb(double _Complex a, int n) {		//[n k_]の計算
	int i,j,l;
	int number=0,hikaku=k[0];
	long double _Complex re_value=1;
	
	for(l=1; l<P; l++) {		//分母k[l]の最大値と、その時の番号lを取得
		if(hikaku < k[l]) {		//各k[l]の大小の比較
			number=l;
			hikaku=k[l];
		}
	}
//	printf("bug:%d  %d\n", number,hikaku);
	
	if(n==0 || n==1 || n==N-1) return 1;
	else {
		for(i=n; i>hikaku; i--) re_value *= sin(i*M_PI/N);	//分子の[n]!の計算
		
		for(i=number+1; i<P; i++) {							//最大値k[i]より番号が大きいk[j]の計算
			for(j=i; j>=2; j--) re_value /= sin(j*M_PI/N);
		}
		for(i=number-1; i>=0; i--) {						//最大値k[i]より番号が小さいk[j]の計算
			for(j=i; j>=2; j--) re_value /= sin(j*M_PI/N);
		}
	}
	return re_value/cpow(sin(M_PI/N),P-1);
}

double _Complex Q_minus(int n, double _Complex q) {
	int i,l;			//カウント用
	long double _Complex re_value = 1;	//戻り値用の変数宣言
	if(n == 0) return 1;
	else {
		for(l=1; l<=n; l++){
			re_value *= (1 - cpow(q,1-N)*cpow(q,n-l));
		}
	}
	return re_value;
}

double _Complex Q_plus(int n, double _Complex q) {
	int i,l;			//カウント用
	long double _Complex re_value = 1;	//戻り値用の変数宣言
	if(n == 0) return 1;
	else {
		for(l=1; l<=n; l++){
			re_value *= (1 - cpow(q,1+N)*cpow(q,n-l));
		}
	}
	return re_value;
}