#include<cstdio>
#include<cstring>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<algorithm>
#include<iostream>
using namespace std;

const double a=5.0,l=30.0,pi=3.1415926535;
double f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,c1,c2,c3,c4,b1,b2,b3,b4;
struct Q{
	double x,y,z,B;//坐标(x,y,z)及磁场 
}s[10000000];

double atanh(double x){
	return 0.5*log((x+1)/(1-x));
}

void set1(double x,double y,double z,double zz){
	f1=sqrt((x+a/2)*(x+a/2)+(y+a/2)*(y+a/2));
	f2=sqrt((x+a/2)*(x+a/2)+(y-a/2)*(y-a/2));
	f3=sqrt((x-a/2)*(x-a/2)+(y+a/2)*(y+a/2));
	f4=sqrt((x-a/2)*(x-a/2)+(y-a/2)*(y-a/2));
	f5=sqrt((x+a/2)*(x+a/2)+(y+a/2)*(y+a/2)+(zz-z)*(zz-z));
	f6=sqrt((x+a/2)*(x+a/2)+(y-a/2)*(y-a/2)+(zz-z)*(zz-z));
	f7=sqrt((x-a/2)*(x-a/2)+(y+a/2)*(y+a/2)+(zz-z)*(zz-z));
	f8=sqrt((x-a/2)*(x-a/2)+(y-a/2)*(y-a/2)+(zz-z)*(zz-z));
	f9=sqrt((x+a/2)*(x+a/2)+(zz-z)*(zz-z));
	f10=sqrt((x-a/2)*(x-a/2)+(zz-z)*(zz-z));
	f11=sqrt((y+a/2)*(y+a/2)+(zz-z)*(zz-z));
	f12=sqrt((y-a/2)*(y-a/2)+(zz-z)*(zz-z));
}

void set2(double x,double y){
	b1=abs(y-a/2);
	b2=abs(y+a/2);
	b3=abs(x-a/2);
	b4=abs(x+a/2);
	c1=sqrt(b4*b4+b1*b1);
	c2=sqrt(b1*b1+b3*b3);
	c3=sqrt(b4*b4+b2*b2);
	c4=sqrt(b2*b2+b3*b3);
}

double ff(int i,int j,double z,double zz){
	double f,b,c;
	if(i==1)c=c1;
	if(i==2)c=c2;
	if(i==3)c=c3;
	if(i==4)c=c4;
	if(j==1)b=b1;
	if(j==2)b=b2;
	if(j==3)b=b3;
	if(j==4)b=b4;
	f=sqrt(1/(c*c-b*b)/(b*b))*atan((zz-z)/sqrt((zz-z)*(zz-z)+c*c)*sqrt((c*c-b*b)/(b*b)));
	return f;
}

double magneticfieldx(double x,double y,double z){ 
	double Bx;
	set1(x,y,z,a/2);
	Bx=(-y-a/2)*atanh((f1-f5)/(abs(x+a/2)-f9))+(-y+a/2)*atanh((f2-f6)/(abs(x+a/2)-f9))+(-y-a/2)*atanh((f3-f7)/(abs(x-a/2)-f10))+(-y+a/2)*atanh((f4-f8)/(abs(x-a/2)-f10));
	set1(x,y,z,-a/2);
	Bx-=(-y-a/2)*atanh((f1-f5)/(abs(x+a/2)-f9))+(-y+a/2)*atanh((f2-f6)/(abs(x+a/2)-f9))+(-y-a/2)*atanh((f3-f7)/(abs(x-a/2)-f10))+(-y+a/2)*atanh((f4-f8)/(abs(x-a/2)-f10));
	return Bx;
}

double magneticfieldy(double x,double y,double z){
	double By;
	set1(x,y,z,a/2);
	By=(x+a/2)*atanh((f2-f6)/(abs(y-a/2)-f12))+(-x+a/2)*atanh((f4-f8)/(abs(y-a/2)-f12))+(x+a/2)*atanh((f1-f5)/(abs(y+a/2)-f11))+(-x+a/2)*atanh((f3-f7)/(abs(y+a/2)-f11));
	set1(x,y,z,-a/2);
	By-=(x+a/2)*atanh((f2-f6)/(abs(y-a/2)-f12))+(-x+a/2)*atanh((f4-f8)/(abs(y-a/2)-f12))+(x+a/2)*atanh((f1-f5)/(abs(y+a/2)-f11))+(-x+a/2)*atanh((f3-f7)/(abs(y+a/2)-f11));
	return By;
}

double magneticfieldz(double x,double y,double z){
	double Bz;
	set2(x,y);
	Bz=(y-a/2)*(x+a/2)*ff(1,1,z,a/2);
	Bz+=(y-a/2)*(-x+a/2)*ff(2,1,z,a/2);
	Bz+=(y+a/2)*(-x-a/2)*ff(3,2,z,a/2);
	Bz+=(y+a/2)*(x-a/2)*ff(4,2,z,a/2);
	Bz+=(y+a/2)*(-x+a/2)*ff(4,3,z,a/2);
	Bz+=(y-a/2)*(x-a/2)*ff(2,3,z,a/2);
	Bz+=(y+a/2)*(-x-a/2)*ff(3,4,z,a/2);
	Bz+=(y-a/2)*(x+a/2)*ff(1,4,z,a/2);
	Bz-=(y-a/2)*(x+a/2)*ff(1,1,z,-a/2);
	Bz-=(y-a/2)*(-x+a/2)*ff(2,1,z,-a/2);
	Bz-=(y+a/2)*(-x-a/2)*ff(3,2,z,-a/2);
	Bz-=(y+a/2)*(x-a/2)*ff(4,2,z,-a/2);
	Bz-=(y+a/2)*(-x+a/2)*ff(4,3,z,-a/2);
	Bz-=(y-a/2)*(x-a/2)*ff(2,3,z,-a/2);
	Bz-=(y+a/2)*(-x-a/2)*ff(3,4,z,-a/2);
	Bz-=(y-a/2)*(x+a/2)*ff(1,4,z,-a/2);
	return Bz*a/2;
}

int main(){
	int i,n;
	double r,x,y,z,Bx,By,Bz,Bx0=0.0,By0=0.0,Bz0=0.0,theta,fai;
	// Bx,By,Bz为在单个磁铁坐标架下的磁场
	//Bx0,By0,Bz0为在磁环坐标架下的磁场
	freopen("a.out","w",stdout);
	x=0.0;y=0;int t=1;
	for(y=-50;y<=50;y+=0.5){
		for(z=-50;z<=50;z+=0.5){
		    Bx0=0;By0=0;Bz0=0;
		    //r=sqrt(y*y+z*z);
		    //if(r>50||(r>27.5&&r<32.4))continue;
			for(i=1;i<=24;i++){//叠加n个磁铁的磁场
			    theta=pi/2+(i-1)*5*pi/12; //即θ角
				fai=(i-1)*pi/12;//即φ角
				Bx=magneticfieldx(x,cos(theta)*(y-l*cos(fai))-sin(theta)*(z-l*sin(fai)),sin(theta)*(y-l*cos(fai))+cos(theta)*(z-l*sin(fai)));
	  		    By=magneticfieldy(x,cos(theta)*(y-l*cos(fai))-sin(theta)*(z-l*sin(fai)),sin(theta)*(y-l*cos(fai))+cos(theta)*(z-l*sin(fai)));
	   	     	Bz=magneticfieldz(x,cos(theta)*(y-l*cos(fai))-sin(theta)*(z-l*sin(fai)),sin(theta)*(y-l*cos(fai))+cos(theta)*(z-l*sin(fai)));
	  		    Bx0+=Bx;
	        	By0+=sqrt(By*By+Bz*Bz)*sin(atan(Bz/By)-theta);
	            Bz0+=sqrt(By*By+Bz*Bz)*sin(atan(Bz/By)-theta);
	        }
			s[t].x=x;
			s[t].y=y;
			s[t].z=z;
			s[t].B=sqrt(Bx0*Bx0+By0*By0+Bz0*Bz0);
			t++;
		}
	}
	for(i=1;i<t;i++){
		printf("%lf\n",s[i].B);
	}printf("\n");cout<<t;
	return 0;
}

