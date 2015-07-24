#include<stdio.h>
#include<complex.h>
#include<math.h>

    /*グローバル変数*/
short p,t,n,i,j,k,l;
short m[30]={0};

/*-------関数のプロトタイプ宣言----------*/
short am(short p_);
short b1m(short p_);
short b2m(void);
double _Complex X(short p_, double _Complex q);
short sigma(short s);
short r(short s);
short r_(short s);
short ik(short s);
double _Complex Q1(short p_,short k,double _Complex q);
double _Complex Q2(short p_,double _Complex q);
short sign(short s);



/*--------メイン関数--------*/
int main(void) {
    short p_;
    double _Complex q;
    long double _Complex J=0, Jall=0;
    printf("Please type value p,t\n");
    scanf("%hd %hd",&p,&t);

    p_ = (p-1)/2;
    for(n=2;n<500;n++){
        Jall=0,J=0;
        //三葉結び目の時の特殊値の計算
        if(p==3){
            for(k=0;k<=p_;k++){ m[k]=0;}
            q = cexp(2*M_PI*I/n);
            while(1){
                J = X(p_,q);    //Two-bridgeによる特殊値の計算
                Jall = Jall + J;
                m[p_]++;
                if(m[p_] > n-1)    break;
            }
            printf("Jall:%.21lf\t",cabs(Jall));
        }

        //それ以外の結び目の特殊値の計算
        else{
            for(k=0;k<=p_;k++){ m[k]=0;}
            q = cexp(2*M_PI*I/n);
            while(1){/*
                printf("%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd\n"
                ,m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8],m[9],m[10],m[11],m[12]
                ,m[13],m[14],m[15],m[16],m[17],m[18],m[19],m[20]);*/
                J = cpow(q,(b1m(p_) + b2m())) * X(p_,q);    //Two-bridgeによる特殊値の計算
                Jall = Jall + J;
                if(m[1] == n-1) break;
                else if(m[p_] < n-1) m[p_]++;
                else if(m[p_] == n-1) {
                  for(j=1; j<=p_-1; j++){
                    if(m[p_-j] < n-1){
                      m[p_-j]++;
                      for(k=p_-j+1; k<=p_; k++){
                        m[k]=m[p_-j];
                      }
                      break;
                    }
                  }
                }
            }
        }
        printf("%hd:Jones:%.32lf\n",n,log(cabs(Jall))/n);
    }
    return 0;
}

/*-----------サブ関数：a(m)の計算----------------*/
short am(short p_){

    short a1=0,a2=0,a3=0,a4=0,A=0;
    for(j=1; j<=p_; j++){
        for(k=r_(j); k<=p_; k++){
            a1=a1+(sigma(ik(k))+sigma(ik(p_+1-k)));
        }
        A = A + a1 * (m[j] - m[j-1]);
        a1=0;
    }
    A=(-1)*A/2;
    
    for(j=1; j<=p_-1; j++){
        a2=a2+(sigma(j+1)+sigma(p_+1-j)) * m[j];
    }

    a3 = ((sigma(p_)+1) / 2) * m[p_];

    for(j=1; j<=p_; j++){
        a4 = a4 + sigma(j);
    }
//    printf("bug-a:%hd\n",A - a2/2 - a3 - a4);
    return A - a2/2 - a3 - a4;
}

    /*サブ関数：b1(m)の計算*/
short b1m(short p_){
    short i,b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,B1=0,B2=0;
    for(k=1; k<=(p-t)/2; k++){
        b1=b1+(1-sigma(ik(k)))/2*(m[ik(k)-1]);
    }
    
    for(k=((p-t)/2)+1;k<=p_; k++){
        b2=b2+(m[ik(k)-1]-(1+sigma(ik(k)))/2*m[ik(k)]);
    }
    
    b3 = (1+sigma(p_))*m[p_];

    for(j=1; j<=p_-1; j++){
        b4=b4+((sigma(j+1)-sigma(j))*m[j]);
    }

    for(k=1; k<=p_-1; k++){
        for(i=k+1; i<=p_; i++){
            b5=b5+(1+sign(ik(k)-ik(i)))/2*(sigma(ik(k))-sigma(ik(i)))*(m[ik(k)]-m[ik(k)-1])*(m[ik(i)]-m[ik(i)-1]);
        }
        B1=B1+b5;
        b5=0;
    }
    B1 = B1 / 2;

    for(j=1; j<=p_; j++){
        for(k=1; k<=r_(j); k++){
            b6=b6+(m[ik(k)]-m[ik(k)-1]);
        }
        B2=B2+sigma(j)*b6*m[j-1];
        b6=0;
    }
//    printf("bug-b1:%hd\n",(-1)*am(p_) + b1 - b2 - b3 + b4/2 - B1 + B2);
    return (-1)*am(p_) + b1 - b2 - b3 + b4/2 - B1 + B2;
}

/*--------サブ関数：b2(m)の計算--------*/
short b2m(void){
    double u = 0;
    if(p<(2*t)){
        for(k=(p-t)/2+1; k<=(t-1)/2; k++){
            u = u + ((1 + sigma(ik(k))) / 2) * m[ik(k) - 1];
        }
    }
    else if(p>(2*t)){
        for(k=(t+1)/2; k<=(p-t)/2; k++){
            u = u + ((1 + sigma(ik(k))) / 2) * m[ik(k) - 1];
        }
        u = u * (-1);
    }
    return u;
}

/*--------サブ関数：X(m)の計算--------*/
double _Complex X(short p_, double _Complex q){
    double _Complex d=0;
    d = pow(-1,m[p_]) * Q1(p_,n-1,q) * Q1(p_,m[p_],q) / (Q1(p_,n-m[p_]-1,q) * Q2(p_,q));
    return d;
}



/*--------各変数の計算--------*/
short sigma(short s){    //σの値の計算
    short u;
    u = r(s) / abs(r(s));
    return u;
}

short ik(short s){        //Iの添え字の判定
    for(l=1; l<=(p-1)/2;l++){    //jは1〜p'まで
        if(s==r_(l)) return l;
    }
}

short r_(short s){        //r'(j)の計算
    short u;
    u = (abs(r(s))+1)/2;
    return u;
}

short r(short s){        //r(j)の計算
    short u;
    u = ((2*s-1)*t)%(2*p);
    if(u<=(-1)*p) u    = u + 2*p;
    else if(p<=u) u = u - 2*p;
    return u;
}

double _Complex Q1(short p_,short k,double _Complex q){        //(q^σ)_kの計算
    double _Complex d=1;
    l=1;
    if(k == 0) return 1;
    else{
        do{
            d = d * (cpow(cpow(q,sigma(p_)),l)-1);
            l++;
        }while(l<=k);
    }
    return d;
}

double _Complex Q2(short p_,double _Complex q){    //(q^-σ)_kの計算
    double _Complex d=1;
    for(j=1; j<=p_; j++){
        if(m[j]-m[j-1] != 0){
            l=1;
            do{
                d = d * (cpow(cpow(q,(-1)*sigma(j)),l)-1);
                l++;
            }while(l<=(m[j]-m[j-1]));
        }
    }
    return d;
}

short sign(short s){    //シグナル計算
    if(s>0) return 1;
    else if(s<0) return -1;
    else return 0;
}
