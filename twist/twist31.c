#include<stdio.h>
#include<complex.h>
#include<math.h>

/*ーーーーグローバル変数ーーーー*/
int N=0;

/*ーーーープロトタイプ宣言ーーーー*/
double _Complex Q_minus(int n, double _Complex q);
double _Complex Q_plus(int n, double _Complex q);

/*
typedef double _Complex Qnotype(double _Complex, int);
double _Complex (*hoge)(double _Complex, int);
Qnotype huga;
hoge = &Q1;
Qnotype qfunc[5] = {&Q2, &Q1, &Q4};

enum {A, B, C, D};

void foo(int Event)
{
    qfunc[Event](3, 4);
}

foo(A);
foo(B);
int printf(char const*, ...)*/

/*ーーーーメイン関数ーーーー*/
int main(void){
    int n;
    for(N=2; N<200; N++){
        long double _Complex J=0;
        double _Complex q = cexp(2*M_PI*I/N);

        for(n=0; n<=N-1; n++){
            J = J + cpow(q,n) * Q_minus(n,q) * Q_plus(n,q);
        }
//        printf("bug=%.25lf  %.25lf\n",creal(q),cimag(q));
        printf("N=%d：%.60lf\n",N,log(cabs(J))/N);
//        printf("test:%.20lf\n",log(4*sqrt(7))/3);
    }
    return 0;
}

/*ーーーーサブ関数ーーーー*/
double _Complex Q_minus(int n, double _Complex q) {
    int i,l;            //カウント用
    long double _Complex re_value = 1;    //戻り値用の変数宣言
    if(n == 0) return 1;
    else {
        for(l=1; l<=n; l++){
            re_value *= (1 - cpow(q,1-N)*cpow(q,n-l));
        }
    }
    return re_value;
}

double _Complex Q_plus(int n, double _Complex q) {
    int i,l;            //カウント用
    long double _Complex re_value = 1;    //戻り値用の変数宣言
    if(n == 0) return 1;
    else {
        for(l=1; l<=n; l++){
            re_value *= (1 - cpow(q,1+N)*cpow(q,n-l));
        }
    }
    return re_value;
}
