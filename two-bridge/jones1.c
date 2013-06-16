#include<stdio.h>
#include<complex.h>
#include<math.h>

/*-------�O���[�o���ϐ�---------*/
int p,t,n;
int m[30]={0};

/*-------�֐��̃v���g�^�C�v�錾----------*/
int am(int p_);
int b1m(int p_);
int b2m(void);
double _Complex X(int p_, double _Complex q);
int sigma(int s);
int r(int s);
int r_(int s);
int ik(int s);
double _Complex Q1(int p_,int k,double _Complex q);
double _Complex Q2(int p_,double _Complex q);
int sign(int s);

/*--------���C���֐�--------*/
int main(void) {
	int p_;
	int i,j,k;
	double _Complex q;
	long double _Complex J=0, Jall=0;
	printf("p,t,n�̒l����͂��Ă��������B");
	scanf("%d %d %d",&p,&t,&n);

	p_ = (p-1)/2;
	q = cexp(2*M_PI*I/n);
	
	/*�O�t���іڂ̎��̓���l�̌v�Z*/
	if(p==3){
		for(k=0;k<=p_;k++){ m[k]=0;}
		while(1){
		printf("%d %d\n",m[0],m[p_]);
		
			J = X(p_,q);			//Two-bridge�ɂ�����l�̌v�Z
			Jall = Jall + J;
		
			m[p_]++;
			if(m[p_] > n-1) break;
		}
		printf("Jall:%.21lf\n",cabs(Jall));
	}

	/*����ȊO�̌��іڂ̓���l�̌v�Z*/
	else{
		for(k=0;k<=p_;k++){ m[k]=0;}
		while(1){
			printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n"
			,m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8],m[9],m[10],m[11],m[12]
					,m[13],m[14],m[15],m[16],m[17],m[18],m[19],m[20]);
			
			J = cpow(q, am(p_)*n + b1m(p_) + b2m()) * X(p_,q);	//Two-bridge�ɂ�����l�̌v�Z
			Jall = Jall + J;
			
			if(m[1] == n-1) break;					//m[1]��n-1(�ő�l)�Ȃ�ΏI��
			if(m[p_] < n-1) m[p_]++;				//�em[p]��n-1(�ő�l)�łȂ��Ȃ�΁A���[�v
			else {									//�em[p]��n-1(�ő�l)�̂Ƃ��A
				for(j=1; j<=p_-1; j++){				//�قȂ���m[i]�̔�r
					if(m[p_-j] < n-1){				//m[p_-j]��n-1(�ő�l)�łȂ��Ȃ��Ȃ�΁A
						m[p_-j]++;					//m[p_-j]�̒l���C���N�������g���A
						for(k=p_-j+1; k<=p_; k++){
							m[k]=m[p_-j];			//m[p_-j+1]�`m[p_]�܂ł̒l��m[p-j]�ɂ���
						}
						break;
					}
				}
			}
		}
		printf("Jall:%.12lf\t",cabs(Jall));
	}
	printf("n=%d:Jones:%.32lf\n\n",n,log(cabs(Jall))/n);
	return 0;
}


/*--------�T�u�֐��Fa(m)�̌v�Z--------*/
int am(int p_){
	int a1=0,a2=0,a3=0,a4=0,A=0;
	int j,k;
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

	a3 = (sigma(p_)+1) * m[p_] / 2;

	for(j=1; j<=p_; j++){
		a4 = a4 + sigma(j);
	}
//	printf("bug-a:%d\n",A - a2/2 - a3 - a4);
	return A - a2/2 - a3 - a4;
}

/*--------�T�u�֐��Fb1(m)�̌v�Z--------*/
int b1m(int p_){
	int b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,B1=0,B2=0;
	int i,j,k;
	for(k=1; k<=(p-t)/2; k++){
		b1=b1+(1-sigma(ik(k)))/2*(m[ik(k)-1]);
	}
//	printf("bug:%d\n",b1);
	
	for(k=((p-t)/2)+1;k<=p_; k++){
		b2 = b2 + (m[ik(k)-1] - (1+sigma(ik(k)))/2*m[ik(k)]);
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

	for(j=1; j<=p_; j++){
		for(k=1; k<=r_(j); k++){
			b6=b6+(m[ik(k)]-m[ik(k)-1]);
		}
		B2=B2+sigma(j)*b6*m[j-1];
		b6=0;
	}
//	printf("bug-b1:%d\n",(-1)*am(p_) + b1 - b2 - b3 + b4/2 - B1 + B2);
	return (-1)*am(p_) + b1 - b2 - b3 + b4/2 - B1/2 + B2;
}

/*--------�T�u�֐��Fb2(m)�̌v�Z--------*/
int b2m(void){
	int k;
	double u = 0;
	if(p<(2*t)){
		for(k=(p-t)/2+1; k<=(t-1)/2; k++){
			u = u + ((1 + sigma(ik(k))) / 2) * m[ik(k) - 1];
		}
	}
	else if(p>(2*t)){
		for(k=(t+1)/2+1; k<=(p-t)/2; k++){
			u = u + ((1 + sigma(ik(k))) / 2) * m[ik(k) - 1];
		}
		u = u * (-1);
	}
//	printf("bug-b2:%.12lf\n",u);
	return u;
}

/*--------�T�u�֐��FX(m)�̌v�Z--------*/
double _Complex X(int p_, double _Complex q){
	double _Complex d=0;
	d = pow(-1,m[p_]) * Q1(p_,n-1,q) * Q1(p_,m[p_],q) / (Q1(p_,n-m[p_]-1,q) * Q2(p_,q));
	return d;
}



/*--------�e�ϐ��̌v�Z--------*/
int sigma(int s){	//�Ђ̒l�̌v�Z
	int u;
	u = r(s) / abs(r(s));
	return u;
}

int ik(int s){		//I�̓Y�����̔���
	int l;
	for(l=1; l<=(p-1)/2;l++){	//j��1�`p'�܂�
		if(s==r_(l)) return l;
	}
}

int r_(int s){		//r'(j)�̌v�Z
	int u;
	u = (abs(r(s))+1)/2;
	return u;
}

int r(int s){		//r(j)�̌v�Z
	int u;
	u = ((2*s-1)*t)%(2*p);
	if(u<=(-1)*p) u	= u + 2*p;
	else if(p<=u) u = u - 2*p;
	return u;
}

double _Complex Q1(int p_,int k,double _Complex q){		//(q^��)_k�̌v�Z
	int l=1;
	double _Complex d=1;
	if(k == 0) return 1;
	else{
		do{
			d = d * (cpow(cpow(q,sigma(p_)),l)-1);
			l++;
		}while(l<=k);
	}
	return d;
}

double _Complex Q2(int p_,double _Complex q){	//(q^-��)_k�̌v�Z
	int j,l;
	double _Complex d=1;
	for(j=1; j<=p_; j++){
		l=1;
		if(m[j]-m[j-1] != 0){
			do{
				d = d * (cpow(cpow(q,(-1)*sigma(j)),l)-1);
				l++;
			}while(l<=(m[j]-m[j-1]));
		}
	}
	return d;
}

int sign(int s){	//�V�O�i���v�Z
	if(s>0) return 1;
	else if(s<0) return -1;
	else return 0;
}
