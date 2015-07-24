/*9_18結び目(41,17)のr(j),r'(j)の計算*/

#include<stdio.h>
#include<math.h>

int main (void) {
    int j,k,l,n,value=0,sum;
    int i[20]={0};
    int sigma[20]={0};
    int sigma_i[20]={0};

    for(j=0; j<=20; j++){
        value = (2*j-1)*17;
        while(1){
            if(value<-41){    value = value + 82;}
            else if(41<value){    value = value - 82;}
            else { break;}
        }
        sigma[j] = value / abs(value);        //σ_jの配列
        i[(abs(value)+1)/2]=j;                //i_r'(j)の配列
//        printf("r(%2d):%3d - r'(%2d):%2d - σ_%2d:%2d\n",j,value,j,(abs(value)+1)/2,j,sigma[j]);
    }
//    printf("\n");

    for(k=0; k<=20; k++){
        sigma_i[k]=sigma[i[k]];
//        printf("σ_i_%2d:%2d\n",k,sigma_i[k]);
    }

    for(n=1; n<=20; n++){
        sum = 0;
        for(l=n; l<=20; l++){
            sum += sigma_i[l] + sigma_i[21-l];
        }
        printf("sum_%2d:%2d\n",n,sum);
    }
    return 0;
}
