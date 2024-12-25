// fft implementation
#include<stdio.h>
#include<math.h> // trigs

void _fft(float buf[],float out[],float exps[],int n,int step){
	// cooley-tukey
	if(step>=n)return;float tr,ti;
	_fft(out,buf,exps,n,step*2);
	_fft(out+step*2,buf+step*2,exps,n,step*2);
	for(int i=0;i<n;i+=2*step){
		tr=exps[2*i]*out[2*(i+step)]-exps[2*i+1]*out[2*(i+step)+1];
		ti=exps[2*i]*out[2*(i+step)+1]+exps[2*i+1]*out[2*(i+step)];
		buf[i]=out[2*i]+tr,buf[i+1]=out[2*i+1]+ti;
		buf[i+n]=out[2*i]-tr,buf[i+n+1]=out[2*i+1]-ti;
	}
}

void _genexps(float buf[],int n,char inv){
	//for(int i=0;i<n;i++)sincosf(M_PI*i*2/n *(inv?-1:1),exps+2*i+1,exps+2*i);
	float a,mlt=inv?-1:1;
	for(int i=0;i<n;i++)a=M_PI*i/n*2,buf[2*i]=cosf(a),buf[2*i+1]=sinf(a)*mlt;
}

void cfft(float buf[],int n,char inv){
	float out[2*n],exps[2*n];
	for(int i=0;i<2*n;i++)out[i]=buf[i];
	_genexps(exps,n,inv); // separate function - allows caching for batch operations
	_fft(buf,out,exps,n,1);
}

void fft(float buf[],int n){
	int s=n;
	if(n&(n-1)){int pow=0,nt=n;
		while(nt)nt>>=1,pow++;s=1<<pow;
	}float in[2*s];
	for(int i=0;i<s;i++)in[2*i]=buf[i%n],in[2*i+1]=0;
	cfft(in,s,0);
	for(int i=0;i<n;i++)buf[i]=sqrt(in[2*i]*in[2*i]+in[2*i+1]*in[2*i+1]);
}
