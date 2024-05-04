#include <stdio.h>
#include <math.h>
#include "lib.h"
#include "mydef.h"
#include "phy_const.h"
#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
#define min(a,b) ((a)<(b)?(a):(b))
//==============constants===============
double dt asm("dt");
const double COE0=q_e/m_e/c_0;
double az asm("az")=1e-100; 				//numeric constants, ~0
double n24 asm("n24")=24;
double ae asm("ae")=1e-12;
double ae2 asm("ae2")=1e-12;
double rcp12 asm("rcp12")=1/12.;
//===========static variables===========
volatile double vlt=1.0,vlt1,vlt0=0.0;				//variables used to ensure the compiler does not ignore some code due to optimization
double t,tau_ asm("tau"),J[4][5],res[4],tau_min,tau_max;	//time, matrix and vector used to solve linear system, and some intermediate variables
double coef0 asm("cf0"),coef1 asm("cf1"),coef2 asm("cf2"),coef3 asm("cf3"),coef4 asm("cf4"),gm;	//intermediate coefficients
Vec3 E asm("e"),B asm("b"),U;				//electric field, magnetic field, 3-momentum
VecN u0 asm("u"),u1 asm("u1"),p asm("p");		//4-velocity in time step n, n+1, 4-position
MatN F,FF,FFF;						//F: electromagnetic tensor 
static Mat3 IB_inv;
//============functions=================
inline void FG(){
	F.x[0][0]=      0;F.x[0][1]= B.x[2];F.x[0][2]=-B.x[1];F.x[0][3]=E.x[0];
	F.x[1][0]=-B.x[2];F.x[1][1]=      0;F.x[1][2]= B.x[0];F.x[1][3]=E.x[1];
	F.x[2][0]= B.x[1];F.x[2][1]=-B.x[0];F.x[2][2]=      0;F.x[2][3]=E.x[2];
	F.x[3][0]= E.x[0];F.x[3][1]= E.x[1];F.x[3][2]= E.x[2];F.x[3][3]=     0;
}
inline void field(Vec3 &E,Vec3 &B,const double t){	//Electromagnetic fields.
	E.x[0]=0; E.x[1]=0; E.x[2]=0; B.x[0]=0; B.x[1]=0; B.x[2]=-3*COE0;
#if(defined(c_txt)||defined(d_txt))
	E.x[2]=COE0;
#endif
	FG();						//calculate electromagnetic tensor F
	E.x[0]=F.x[0][3]; E.x[1]=F.x[1][3]; E.x[2]=F.x[2][3]; B.x[0]=F.x[1][2]; B.x[1]=-F.x[0][2]; B.x[2]=F.x[0][1];	//field in computational frame
}
inline void Comp_coe(){					//compute coefficients
	double BB=B*B,BE=B*E,BU=B*U;coef0=2.*(BE)*(BU)+(BE)*(BE)*dt,coef1=2.*(U*(B.CPM()*E)+BB*gm),coef2=2.*E*U+E*E*dt-BB*dt,coef3=2.*gm,coef4=-dt;
}
inline double Q_f(double tau){				//quartic function
	return (((coef0*tau+coef1)*tau+coef2)*tau+coef3)*tau+coef4;
}
inline double D_qf(double &tau){				//derivative of the quartic function
	return ((4.*coef0*tau+3.*coef1)*tau+2.*coef2)*tau+coef3;
}
inline void Tau_Max(){										//calculate the upper limit of \tau
	tau_min=0;
	double BB=B*B,EE=E*E,BE=B*E,TEMP=EE-BB,Res;
	if(fabs(TEMP)*1e-12<fabs(BE)) Res=sqrt(2/(TEMP+sqrt(TEMP*TEMP+4*BE*BE)));		//when BE*BE is nonzero
	else if(EE>BB) Res=sqrt(1./(TEMP));							//when BE*BE is approximately zero and EE>BB
	else Res=dt/(1+gm);
	tau_max=min(Res,dt/1+gm);
}
inline void Ens_conv(){					//ensure that f is a convex function in [\tau_{min},\tau_{max}]
	if(fabs(coef0)>1e-100){
		double temp=9.*coef1*coef1-24.*coef0*coef2;
		if(temp>=0){
			temp=sqrt(temp);
			double t1=(-3.0*coef1+temp)/(12.*coef0),t2=(-3.0*coef1-temp)/(12.*coef0);
			if(tau_min<t1&&t1<tau_max){if(Q_f(t1)>0) tau_max=t1; else tau_min=t1;}
			if(tau_min<t2&&t2<tau_max){if(Q_f(t2)>0) tau_max=t2; else tau_min=t2;}
		}
	}
	else if((fabs(coef1)>1e-100)){
		double t1=-coef2/(3.*coef1);
		if(tau_min<t1&&t1<tau_max){if(Q_f(t1)>0) tau_max=t1; else tau_min=t1;}
	}
}
inline void Tau_Newt(){		//Newton iteration to solve \tau
	double a=D_qf(tau_min),b=D_qf(tau_max),dt;
	if(a>b){dt=Q_f(tau_min)/a;tau_=tau_min-dt;} else{dt=Q_f(tau_max)/b;tau_=tau_max-dt;}
	while(likely(fabs(dt)>1e-12*tau_)){dt=Q_f(tau_)/D_qf(tau_);tau_-=dt;}
}
inline void Comp_J(){		//compute augmented Jacobian matrix
	J[0][0]=1;J[0][1]=-B.x[2]*tau_;J[0][2]=B.x[1]*tau_;J[0][3]=-E.x[0]*tau_;J[0][4]=u0.x[0]-B.x[1]*u0.x[2]*tau_+B.x[2]*u0.x[1]*tau_+E.x[0]*u0.x[3]*tau_;
	J[1][0]=B.x[2]*tau_;J[1][1]=1;J[1][2]=-B.x[0]*tau_;J[1][3]=-E.x[1]*tau_;J[1][4]=u0.x[1]-B.x[2]*u0.x[0]*tau_+B.x[0]*u0.x[2]*tau_+E.x[1]*u0.x[3]*tau_;
	J[2][0]=-B.x[1]*tau_;J[2][1]=B.x[0]*tau_;J[2][2]=1;J[2][3]=-E.x[2]*tau_;J[2][4]=u0.x[2]-B.x[0]*u0.x[1]*tau_+B.x[1]*u0.x[0]*tau_+E.x[2]*u0.x[3]*tau_;
	J[3][0]=-E.x[0]*tau_;J[3][1]=-E.x[1]*tau_;J[3][2]=-E.x[2]*tau_;J[3][3]=1;  J[3][4]=u0.x[3]+E.x[0]*u0.x[0]*tau_+E.x[1]*u0.x[1]*tau_+E.x[2]*u0.x[2]*tau_;
}
inline void Gauss_Jordan(){	//Gauss Jordan elimination with complete pivoting
	const int N=4;
	int J_order[N];
	for(int i=0;i<N;i++) J_order[i]=i;
	for(int k=0;k<N;k++) {
		int tag_i,tag_j,m;
		double temp=0.;
		for(int i=k;i<N;i++) for(int j=k;j<N;j++){if(fabs(J[i][j])>temp){temp=fabs(J[i][j]);tag_i=i;tag_j=j;}}
		if(unlikely(temp<1e-30)){fprintf(stderr,"NOT full rank, rank= %d, order %d\n",k,J_order[k]);break;}
		m=J_order[k];J_order[k]=J_order[tag_j];J_order[tag_j]=m;
		for(int i=k;i<N+1;i++){temp=J[k][i];J[k][i]=J[tag_i][i];J[tag_i][i]=temp;}
		for(int i=0;i<N;i++){temp=J[i][k];J[i][k]=J[i][tag_j];J[i][tag_j]=temp;}
		temp=J[k][k];
		for(int i=k;i<N+1;i++)J[k][i]/=temp;
			for(int i=0;i<N;i++) if(i!=k){for(int j=k+1;j<N+1;j++) J[i][j]-=J[k][j]*J[i][k]; J[i][k]=0;}	
	}
	for(int i=0;i<N;i++) res[J_order[i]]=J[i][N];
}
void Comp_tau(){
	field(E,B,t+0.5*dt);				//compute field
	Comp_coe();
	Tau_Max();
	Ens_conv();
	Tau_Newt();
}
inline void Update_imp1(){				//algorithm(implementation 1)	
	Comp_J();
	Gauss_Jordan();
	u1=u0;
	u0.x[0]=U.x[0]=res[0];u0.x[1]=U.x[1]=res[1];u0.x[2]=U.x[2]=res[2];u0.x[3]=gm=res[3];
	p+=(u0+u1)*tau_;
}
inline void Comp_ib_inv(){
	double t2=tau_*tau_;
	IB_inv.x[0][0]=1.+B.x[0]*B.x[0]*t2;		IB_inv.x[0][1]=B.x[2]*tau_+B.x[0]*B.x[1]*t2;	IB_inv.x[0][2]= -B.x[1]*tau_+B.x[0]*B.x[2]*t2;
	IB_inv.x[1][0]= -B.x[2]*tau_+B.x[0]*B.x[1]*t2;	IB_inv.x[1][1]=1+B.x[1]*B.x[1]*t2;		IB_inv.x[1][2]=B.x[0]*tau_+B.x[1]*B.x[2]*t2;
	IB_inv.x[2][0]=B.x[1]*tau_+B.x[0]*B.x[2]*t2;	IB_inv.x[2][1]= -B.x[0]*tau_+B.x[1]*B.x[2]*t2;	IB_inv.x[2][2]=1.+B.x[2]*B.x[2]*t2;
	IB_inv=IB_inv*(1./(1.+B*B*t2));
}
inline void Update_imp2(){
	u1=u0;
	Comp_ib_inv();
	U=IB_inv*(2.*U+E*dt)-U;
	gm=sqrt(1+U*U);u0.x[0]=U.x[0];u0.x[1]=U.x[1];u0.x[2]=U.x[2];u0.x[3]=gm;
	p+=(u0+u1)*tau_;				//update 4-position
}

void Asm_opt(){
	field(E,B,t+0.5*dt);									//compute field
	asm(
		"vmovsd b(%rip),%xmm0\n"		//B0
		"vmovsd 8+b(%rip),%xmm1\n"		//B1
		"vmovsd 16+b(%rip),%xmm2\n"		//B2
		"vmovsd e(%rip),%xmm3\n"		//E0
		"vmovsd 8+e(%rip),%xmm4\n"		//E1
		"vmovsd 16+e(%rip),%xmm5\n"		//E2
		"vmovsd u(%rip),%xmm6\n"		//U0
		"vmovsd 8+u(%rip),%xmm7\n"		//U1
		"vmovsd 16+u(%rip),%xmm8\n"		//U2
		"vmovsd 24+u(%rip),%xmm15\n"		//gm
		"vmovsd dt(%rip),%xmm14\n"		//dt
		"vpcmpeqb %ymm9,%ymm9,%ymm9\n"		//ymm9=(1,1,1,1)
		"vpsrlq $0x36,%ymm9,%ymm9\n"		//ymm9=(1,1,1,1)
		"vpsllq $0x34,%ymm9,%ymm9\n"		//ymm9=(1,1,1,1)
		"vaddsd %xmm9,%xmm15,%xmm9\n"		//ymm9=(1+gm)
		"vdivsd %xmm9,%xmm14,%xmm13\n"		//xmm13=TEMP2=dt/(1+gm)
		"vmulsd %xmm2,%xmm4,%xmm9\n"		//BXE
		"vfmsub231sd %xmm1,%xmm5,%xmm9\n"	//BXE
		"vmulsd %xmm0,%xmm5,%xmm10\n"		//BXE
		"vfmsub231sd %xmm2,%xmm3,%xmm10\n"	//BXE
		"vmulsd %xmm1,%xmm3,%xmm11\n"		//BXE
		"vfmsub231sd %xmm0,%xmm4,%xmm11\n"	//BXE
		"vmulsd %xmm9,%xmm6,%xmm12\n"		//U(BXE)
		"vfmadd231sd %xmm10,%xmm7,%xmm12\n"	//U(BXE)
		"vfmadd231sd %xmm11,%xmm8,%xmm12\n"	//xmm12=U(BXE)
		"vmulsd %xmm0,%xmm6,%xmm11\n"		//BU
		"vfmadd231sd %xmm1,%xmm7,%xmm11\n"	//BU
		"vfmadd231sd %xmm2,%xmm8,%xmm11\n"	//xmm11=BU
		"vmulsd %xmm0,%xmm0,%xmm9\n"		//BB
		"vfmadd231sd %xmm1,%xmm1,%xmm9\n"	//BB
		"vfmadd231sd %xmm2,%xmm2,%xmm9\n"	//xmm9=BB
		"vmulsd %xmm0,%xmm3,%xmm10\n"		//BE
		"vmovsd dt(%rip),%xmm0\n"		//dt
		"vfmadd231sd %xmm1,%xmm4,%xmm10\n"	//BE
		"vfmadd231sd %xmm2,%xmm5,%xmm10\n"	//xmm10=BE
		"vmulsd %xmm5,%xmm8,%xmm8\n"		//EU
		"vfmadd231sd %xmm3,%xmm6,%xmm8\n"	//EU
		"vfmadd231sd %xmm4,%xmm7,%xmm8\n"	//xmm8=EU
		"vaddsd %xmm8,%xmm8,%xmm2\n"		//coef2
		"vmulsd %xmm3,%xmm3,%xmm7\n"		//EE
		"vfmadd231sd %xmm4,%xmm4,%xmm7\n"	//EE
		"vfmadd231sd %xmm5,%xmm5,%xmm7\n"	//xmm7=EE
		"vsubsd %xmm9,%xmm7,%xmm8\n"		//xmm8=TEMP0=EE-BB
		"vmovsd ae2(%rip),%xmm7\n"		//xmm7=1e-12
		"vmulsd %xmm10,%xmm10,%xmm6\n"		//xmm6=TEMP1=BE*BE
		"vxorpd %xmm4,%xmm4,%xmm4\n"		//coef4
		"vsubsd %xmm14,%xmm4,%xmm4\n"		//coef4=-dt;
		"vfmadd231sd %xmm14,%xmm8,%xmm2\n"	//xmm2=coef2=2.*E*U+EE*dt-BB*dt
		"vaddsd %xmm15,%xmm15,%xmm3\n"		//xmm3=coef3=2*gm
		"vmulsd %xmm10,%xmm11,%xmm0\n"		//coef0
		"vaddsd %xmm0,%xmm0,%xmm0\n"		//coef0
		"vfnmadd231sd %xmm4,%xmm6,%xmm0\n"	//xmm0=coef0
		"vfmadd231sd %xmm15,%xmm9,%xmm12\n"	//coef1
		"vpcmpeqb %ymm9,%ymm9,%ymm9\n"		//ymm9=(1,1,1,1)
		"vpsrlq $0x36,%ymm9,%ymm9\n"		//ymm9=(1,1,1,1)
		"vpsllq $0x34,%ymm9,%ymm9\n"		//ymm9=(1,1,1,1)
		"vaddsd %xmm12,%xmm12,%xmm1\n"		//xmm1=coef1
		"vmovsd az(%rip),%xmm12\n"		//xmm12=1e-100
		"vmovapd %ymm13,%ymm11\n"		//xmm11=TEMP3=TEMP2
		"vmulsd %xmm8,%xmm7,%xmm14\n"		//xmm14=TEMP0*1e-12;
		"vpsllq $0x1,%ymm14,%ymm14\n"		//xmm14=fabs(EE-BB)*1e-12
		"vpsrlq $0x1,%ymm14,%ymm14\n"		//xmm14=fabs(EE-BB)*1e-12
		"vmovapd %ymm10,%ymm15\n"		//ymm15=BE
		"vpsllq $0x1,%ymm15,%ymm15\n"		//ymm15=fabs(BE)
		"vpsrlq $0x1,%ymm15,%ymm15\n"		//ymm15=fabs(BE)
		"vcomisd %xmm15,%xmm14\n"		//if(fabs(EE-BB)*1e-12<fabs(BE))	//when BE*BE is nonzero
		"jae .dst1\n"
		"vaddsd %xmm6,%xmm6,%xmm11\n"		//2*TEMP1  TEMP3=sqrt(2./(TEMP0+sqrt(TEMP0*TEMP0+4*TEMP1)));
		"vaddsd %xmm11,%xmm11,%xmm11\n"		//4*TEMP1
		"vfmadd231sd %xmm8,%xmm8,%xmm11\n"	///xmm11=TEMP0*TEMP0+4*TEMP1;
		"vsqrtpd %xmm11,%xmm11\n"		//sqrt(xmm11)
		"vaddsd %xmm8,%xmm11,%xmm11\n"		//TEMP0+sqrt(TEMP0*TEMP0+4*TEMP1);
		"vaddsd %xmm9,%xmm9,%xmm7\n"		//ymm7=(2,2,2,2)
		"vdivsd %xmm11,%xmm7,%xmm11\n"		//xmm11=2./(TEMP0+sqrt(TEMP0*TEMP0+4*TEMP1))
		"vsqrtpd %xmm11,%xmm11\n"		//xmm11=TEMP3=sqrt(2./(TEMP0+sqrt(TEMP0*TEMP0+4*TEMP1)));
		"jmp .dst0\n"
		".dst1:\n"
		"vxorpd %xmm7,%xmm7,%xmm7\n"		//0
		"vcomisd %xmm7,%xmm8\n"			//else if(TEMP0>0)	//When BE*BE is approximately zero and EE>BB
		"jbe .dst0\n"
		"vdivsd %xmm8,%xmm9,%xmm11\n"		//TEMP3=sqrt(1./TEMP0);
		"vsqrtpd %xmm11,%xmm11\n"
		".dst0:\n"
		"vminsd %xmm13,%xmm11,%xmm10\n"		//find out tau_max(the upper limit of \tau), xmm10=tau_max=min(res,dt/1+gm)
		"vmovsd n24(%rip),%xmm14\n"		//xmm14=24
		"vmovsd rcp12(%rip),%xmm15\n"		//xmm15=1/12
		"vaddsd %xmm0,%xmm0,%xmm5\n"
		"vaddsd %xmm5,%xmm5,%xmm5\n"		//xmm5=4*coef0
		"vaddsd %xmm1,%xmm1,%xmm6\n"
		"vaddsd %xmm6,%xmm1,%xmm6\n"		//xmm6=3*coef1
		"vaddsd %xmm2,%xmm2,%xmm7\n"		//xmm7=2*coef2
		"vxorpd %xmm9,%xmm9,%xmm9\n"		//tau_min=0
		"vpcmpeqb %ymm11,%ymm11,%ymm11\n"
		"vpsllq $0x3f,%ymm11,%ymm13\n"		//-
		"vpsrlq $0x1,%ymm11,%ymm11\n"		// fabs
		"vandpd %ymm11,%ymm0,%ymm8\n"		//fabs(coef0)
		"comisd %xmm12,%xmm8\n"			//if(fabs(coef0)>1e-100)
		"jbe .dst2\n"
		"vmulsd %xmm14,%xmm0,%xmm8\n"
		"vmulsd %xmm8,%xmm2,%xmm8\n"
		"vfmsub231sd %xmm6,%xmm6,%xmm8\n"	//double temp=9.*coef1*coef1-24.*coef0*coef2;
		"vxorpd %ymm12,%ymm12,%ymm12\n"
		"comisd %xmm12,%xmm8\n"			//if(temp>=0)
		"jb .dst2\n"
		"vsqrtpd %xmm8,%xmm8\n"
		"vdivsd %xmm0,%xmm15,%xmm12\n"		//xmm12=TEMP0=1/(12*coef0);
		"vsubsd %xmm6,%xmm8,%xmm14\n"		//xmm14=t1  double t1=(-3.0*coef1+temp)*TEMP0,t2=(-3.0*coef1-temp)*TEMP0;
		"vmulsd %xmm14,%xmm12,%xmm14\n"
		"vxorpd %xmm12,%xmm13,%xmm12\n"		//xmm12=-1/(12*coef0);
		"vaddsd %xmm6,%xmm8,%xmm15\n"
		"vmulsd %xmm15,%xmm12,%xmm15\n"		//xmm15=t2
		"comisd %xmm14,%xmm9\n"			//if(tau_min<t1&&t1<tau_max){if((((coef0*t1+coef1)*t1+coef2)*t1+coef3)*t1+coef4>0) tau_max=t1; else tau_min=t1;}
		"jae .dst8\n"
		"comisd %xmm10,%xmm14\n"
		"jae .dst8\n"
		"vmovapd %ymm0,%ymm8\n"
		"vfmadd213sd %xmm1,%xmm14,%xmm8\n"
		"vfmadd213sd %xmm2,%xmm14,%xmm8\n"
		"vfmadd213sd %xmm3,%xmm14,%xmm8\n"
		"vfmadd213sd %xmm4,%xmm14,%xmm8\n"	//xmm8=(((coef0*t1+coef1)*t1+coef2)*t1+coef3)*t1+coef4
		"vxorpd %xmm12,%xmm12,%xmm12\n"	
		"comisd %xmm12,%xmm8\n"
		"jbe .dst7\n"
		"vmovapd %ymm14,%ymm10\n"		//if((((coef0*t1+coef1)*t1+coef2)*t1+coef3)*t1+coef4>0) tau_max=t1;
		"jmp .dst8\n"
		".dst7:\n"
		"vmovapd %ymm14,%ymm9\n"		//else tau_min=t1;
		".dst8:\n"
		"comisd %xmm15,%xmm9\n"			//if(tau_min<t2&&t2<tau_max){if((((coef0*t2+coef1)*t2+coef2)*t2+coef3)*t2+coef4>0) tau_max=t2; else tau_min=t2;}
		"jae .dst2\n"
		"comisd %xmm10,%xmm15\n"
		"jae .dst2\n"
		"vmovapd %ymm0,%ymm8\n"
		"vfmadd213sd %xmm1,%xmm15,%xmm8\n"
		"vfmadd213sd %xmm2,%xmm15,%xmm8\n"
		"vfmadd213sd %xmm3,%xmm15,%xmm8\n"
		"vfmadd213sd %xmm4,%xmm15,%xmm8\n"	//xmm8=(((coef0*t2+coef1)*t2+coef2)*t2+coef3)*t2+coef4
		"vxorpd %xmm12,%xmm12,%xmm12\n"	
		"comisd %xmm12,%xmm8\n"
		"jbe .dst9\n"
		"vmovapd %ymm15,%ymm10\n"		//if((((coef0*t2+coef1)*t2+coef2)*t2+coef3)*t2+coef4>0) tau_max=t2;
		"jmp .dst2\n"
		".dst9:\n"
		"vmovapd %ymm15,%ymm9\n"		//else tau_min=t2;
		".dst2:\n"
		"vandpd %ymm11,%ymm1,%ymm8\n"		//fabs(coef1)
		"comisd %xmm12,%xmm8\n"			//if((fabs(coef1)>1e-100))
		"jbe .dst5\n"
		"vdivsd %xmm6,%xmm2,%xmm14\n"		//t1=-coef2/(3.*coef1);
		"vxorpd %xmm14,%xmm13,%xmm14\n"		//t1=-coef2/(3.*coef1);
		"comisd %xmm14,%xmm9\n"			//if(tau_min<t1&&t1<tau_max)
		"jae .dst5\n"
		"comisd %xmm10,%xmm14\n"
		"jae .dst5\n"
		"vmovapd %ymm0,%ymm8\n"
		"vfmadd213sd %xmm1,%xmm14,%xmm8\n"
		"vfmadd213sd %xmm2,%xmm14,%xmm8\n"
		"vfmadd213sd %xmm3,%xmm14,%xmm8\n"
		"vfmadd213sd %xmm4,%xmm14,%xmm8\n"	//xmm8=(((coef0*t1+coef1)*t1+coef2)*t1+coef3)*t1+coef4
		"vxorpd %xmm12,%xmm12,%xmm12\n"
		"comisd %xmm12,%xmm8\n"			//if(Q_f(t1)>0) tau_max=t1; 
		"jbe .dst10\n"
		"vmovapd %ymm14,%ymm10\n"
		"jmp .dst5\n"
		".dst10:\n"
		"vmovapd %ymm14,%ymm9\n"		//else tau_min=t1;	
		".dst5:\n"
		"vmovsd ae(%rip),%xmm13\n"
		"vmovapd %ymm5,%ymm8\n"			//xmm8=a=D_qf(tau_min)
		"vfmadd213sd %xmm6,%xmm9,%xmm8\n"
		"vfmadd213sd %xmm7,%xmm9,%xmm8\n"
		"vfmadd213sd %xmm3,%xmm9,%xmm8\n"
		"vmovapd %ymm5,%ymm12\n"		//xmm12=b=D_qf(tau_max)
		"vfmadd213sd %xmm6,%xmm10,%xmm12\n"
		"vfmadd213sd %xmm7,%xmm10,%xmm12\n"
		"vfmadd213sd %xmm3,%xmm10,%xmm12\n"
		"comisd %xmm12,%xmm8\n"			//if(a>b)
		"jbe .dst11\n"
		"vmovapd %ymm0,%ymm11\n"		//xmm11=Q_f(tau_min)
		"vfmadd213sd %xmm1,%xmm9,%xmm11\n"
		"vfmadd213sd %xmm2,%xmm9,%xmm11\n"
		"vfmadd213sd %xmm3,%xmm9,%xmm11\n"
		"vfmadd213sd %xmm4,%xmm9,%xmm11\n"
		"vdivsd %xmm8,%xmm11,%xmm12\n"		//xmm12=dtau=Q_f(tau_min)/a
		"vsubsd %xmm12,%xmm9,%xmm11\n"		//xmm11=tau-tau_min-dtau;
		"jmp .dst12\n"
		".dst11:\n"				//else{dtau=Q_f(tau_max)/b;tau_=tau_max-dtau;}
		"vmovapd %ymm0,%ymm11\n"		//xmm11=dtau=Q_f(tau_max)
		"vfmadd213sd %xmm1,%xmm10,%xmm11\n"
		"vfmadd213sd %xmm2,%xmm10,%xmm11\n"
		"vfmadd213sd %xmm3,%xmm10,%xmm11\n"
		"vfmadd213sd %xmm4,%xmm10,%xmm11\n"
		"vdivsd %xmm12,%xmm11,%xmm12\n"		//xmm12=dtau=Q_f(tau_max)/b
		"vsubsd %xmm12,%xmm10,%xmm11\n"		//xmm11=tau_=tau_max-dtau;
		".dst12:\n"
		"vpsllq $0x1,%ymm12,%ymm12\n"
		"vpsrlq $0x1,%ymm12,%ymm12\n"		//fabs(dtau)
		"vmulsd %xmm13,%xmm11,%xmm9\n"		//xmm9=1e-12*tau
		"comisd %xmm9,%xmm12\n"	//while(fabs(dtau)>1e-12*tau_){dtau=Q_f(tau_)/D_qf(tau_);tau_-=dtau;}
		"jbe .dst13\n"
		"vmovapd %ymm0,%ymm9\n"
		"vfmadd213sd %xmm1,%xmm11,%xmm9\n"
		"vfmadd213sd %xmm2,%xmm11,%xmm9\n"
		"vfmadd213sd %xmm3,%xmm11,%xmm9\n"
		"vfmadd213sd %xmm4,%xmm11,%xmm9\n"
		"vmovapd %ymm5,%ymm12\n"
		"vfmadd213sd %xmm6,%xmm11,%xmm12\n"
		"vfmadd213sd %xmm7,%xmm11,%xmm12\n"
		"vfmadd213sd %xmm3,%xmm11,%xmm12\n"
		"vdivsd %xmm12,%xmm9,%xmm12\n"
		"vsubsd %xmm12,%xmm11,%xmm11\n"
		"jmp .dst12\n"
		".dst13:\n"
		"vmovsd %xmm11,tau(%rip)\n"
		//update
		"vbroadcastsd %xmm11,%ymm11\n"		//ymm11=(tau,tau,tau,tau)
		"vmulpd b(%rip),%ymm11,%ymm15\n"	//ymm15=B*tau
		"vbroadcastsd dt(%rip),%ymm8\n"		//ymm8=(dt,dt,dt,dt)
		"vmulpd e(%rip),%ymm8,%ymm14\n"		//ymm14=E*dt
		"vmovapd p(%rip),%ymm13\n"		//ymm13=p
		"vmovapd u(%rip),%ymm12\n"		//ymm12=u
		"vpcmpeqb %ymm10,%ymm10,%ymm10\n"	//ymm10=-
		"vpsllq $0x3f,%ymm10,%ymm10\n"		//ymm10=-
		"vpcmpeqb %ymm1,%ymm1,%ymm1\n"		//%%ymm1=(1,1,1,1)
		"vpsrlq $0x36,%ymm1,%ymm1\n"		//%%ymm1=(1,1,1,1)
		"vpsllq $0x34,%ymm1,%ymm1\n"		//%%ymm1=(1,1,1,1)
		"vmovapd %ymm1,%ymm9\n"			//ymm9=1
		"vfmadd231pd %ymm15,%ymm15,%ymm9\n"	//ymm9=1+B0*B0
		"vbroadcastsd %xmm15,%ymm5\n"		//ymm5=B0,B0,B0,B0
		"vextractf128 $0b1,%ymm15,%xmm2\n"	//ymm15=B0,B1,B2,   xmm2=B2,,,
		"vbroadcastsd %xmm2,%ymm7\n"		//ymm7=B2,B2,B2,B2
		"vfmadd231pd %ymm2,%ymm2,%ymm9\n"	//ymm9=1+B0*B0+B2*B2
		"vfmadd231pd %ymm12,%ymm11,%ymm13\n"	//p+=u1*tau
		"vunpckhpd %ymm15,%ymm15,%ymm2\n"	//ymm2=B1,x
		"vbroadcastsd %xmm2,%ymm6\n"		//ymm6 B1,B1,B1,B1
		"vfmadd231pd %ymm6,%ymm6,%ymm9\n"	//ymm9=1+B0*B0+B2*B2+B1*B1
		"vdivpd %ymm9,%ymm1,%ymm9\n"		//ymm9=1./ymm9  coef
		"vaddpd %ymm12,%ymm12,%ymm8\n"		//ymm8=2*U
		"vaddpd %ymm8,%ymm14,%ymm8\n"		//ymm8=2*U+E*dt
		"vbroadcastsd %xmm9,%ymm9\n"		//ymm9=ymm9,ymm9,ymm9,ymm9
		"vxorpd %ymm7,%ymm10,%ymm2\n"		//ymm2=-B2,-B2,-B2,-B2
		"vunpcklpd %ymm2,%ymm1,%ymm2\n"		//xmm2=1,-B2
		"vinsertf128 $0b1,%xmm6,%ymm2,%ymm2\n"	//ymm2=1,-B2,B1,x
		"vfmadd231pd %ymm15,%ymm5,%ymm2\n"	//Mat+=BB dyad
		"vunpcklpd %ymm1,%ymm7,%ymm3\n"		//ymm3=B2,1
		"vxorpd %ymm5,%ymm10,%ymm4\n"		//ymm4=-B0,-B0,-B0,-B0
		"vmulpd %ymm2,%ymm9,%ymm2\n"		//mult coef
		"vinsertf128 $0b1,%xmm4,%ymm3,%ymm3\n"	//ymm3=B2,1,-B0,x
		"vfmadd231pd %ymm15,%ymm6,%ymm3\n"	//Mat+=BB dyad
		"vmulpd %ymm3,%ymm9,%ymm3\n"		//mult coef
		"vxorpd %ymm6,%ymm10,%ymm4\n"		//ymm4=-B1,-B1,-B1,-B1
		"vunpcklpd %ymm5,%ymm4,%ymm4\n"		//ymm4=-B1,B0
		"vbroadcastsd %xmm8,%ymm5\n"		//ymm5 2*U+E*dt comp
		"vunpckhpd %ymm8,%ymm8,%ymm6\n"		//ymm6 2*U+E*dt comp
		"vbroadcastsd %xmm6,%ymm6\n"		//ymm6 2*U+E*dt comp
		"vextractf128 $0b1,%ymm8,%xmm8\n"	//ymm8 2*U+E*dt comp
		"vbroadcastsd %xmm8,%ymm8\n"		//ymm8 2*U+E*dt comp
		"vinsertf128 $0b1,%xmm1,%ymm4,%ymm4\n"	//ymm4=-B1,B0,1,x
		"vfmadd231pd %ymm15,%ymm7,%ymm4\n"	//Mat+=BB dyad
		"vmulpd %ymm4,%ymm9,%ymm4\n"		//mult coef
		"vfmsub231pd %ymm2,%ymm5,%ymm12\n"	//u=IB_inv*(2*U+E*dt)-U
		"vfmadd231pd %ymm3,%ymm6,%ymm12\n"	//u=IB_inv*(2*U+E*dt)-U
		"vfmadd231pd %ymm4,%ymm8,%ymm12\n"	//u=IB_inv*(2*U+E*dt)-U
		"vmovapd %ymm1,%ymm2\n"			//ymm2=1,1,1,1
		"vfmadd231pd %ymm12,%ymm12,%ymm2\n"	//ymm2= , ,1+u2*u2,
		"vperm2f128 $0b0,%ymm12,%ymm12,%ymm5\n"	//ymm5=u0,u1,u0,u1
		"vfmadd231pd %ymm5,%ymm5,%ymm2\n"	//ymm2= , ,1+u2*u2+u0*u0,
		"vpermilpd $0b0,%ymm2,%ymm2\n"		//ymm2= , , ,1+u2*u2+u0*u0
		"vfmadd231pd %ymm5,%ymm5,%ymm2\n"	//ymm2= , , ,1+u2*u2+u0*u0+u1*u1
		"vsqrtpd %ymm2,%ymm2\n"			//ymm2= , , ,gm
		"vblendpd $0b1000,%ymm2,%ymm12,%ymm12\n"//ymm12=u0
		"vfmadd231pd %ymm12,%ymm11,%ymm13\n"	//p+=u1*tau
		"vmovapd %ymm12,u1(%rip)\n"		//save u to u1, u1 is in memory 
		"vmovapd %ymm12,u(%rip)\n"		//save u to u0, u0 is in memory 
		"vmovapd %ymm13,p(%rip)\n"		//save p to memory
	);
}
const double eps_th=1e-3;
double sinc(double x){
	if(likely(abs(x)>eps_th)) return sin(x)/x;
	else return 1.-1./6.*x*x+1./120*x*x*x*x;
}
double sinci(double x){
	if(likely(abs(x)>eps_th)) return sinh(x)/x;	//sinh(x)=(e^x-e^{-x})/2
	else return 1.+1./6.*x*x+1./120*x*x*x*x;	//sinh(x)=conh(x)'=sinh(x)"
}
void anal_solu(){
	field(E,B,t+0.5*dt);				//compute fields
	double TAU,I1=E*E-B*B,I2=E*B;
	double kappa=1/sqrt(2)*sqrt(I1+sqrt(I1*I1+4*I2*I2)),omega=1/sqrt(2)*sqrt(-I1+sqrt(I1*I1+4*I2*I2));
	int tar;
	if(unlikely(kappa<eps_th&&omega<eps_th)){kappa=eps_th;omega=eps_th;tar=1;} else tar=0;
	FF=F*F;
	const MatN I=1.,P_k=(omega*omega*I+FF)/(kappa*kappa+omega*omega),P_o=(kappa*kappa*I-FF)/(kappa*kappa+omega*omega);
	const VecN U_k=P_k*u0,U_o=P_o*u0;
	double tau0=0.,tau1=dt,f0,FUk=F.x[3][0]*U_k.x[0]+F.x[3][1]*U_k.x[1]+F.x[3][2]*U_k.x[2],FUo=F.x[3][0]*U_o.x[0]+F.x[3][1]*U_o.x[1]+F.x[3][2]*U_o.x[2];
	f0=(U_k.x[3]*sinci(kappa*tau0)+U_o.x[3]*sinc(omega*tau0)+.5*FUk*sinci(.5*kappa*tau0)*sinci(.5*kappa*tau0)*tau0+.5*FUo*sinc(.5*omega*tau0)*sinc(.5*omega*tau0)*tau0)*tau0-dt;
/*
	for(;;){					//secant method
		double f1=(U_k.x[3]*sinci(kappa*tau1)+U_o.x[3]*sinc(omega*tau1)+.5*FUk*sinci(.5*kappa*tau1)*sinci(.5*kappa*tau1)*tau1+.5*FUo*sinc(.5*omega*tau1)*sinc(.5*omega*tau1)*tau1)*tau1-dt;
		TAU=tau0-f0/(f1-f0)*(tau1-tau0);
		if(fabs(tau1-tau0)<1e-6*TAU) break;
		f0=f1;tau0=tau1;tau1=TAU;
	}
*/
	for(;;){					//bisection method
		TAU=.5*(tau1+tau0);
		double f1=(U_k.x[3]*sinci(kappa*TAU)+U_o.x[3]*sinc(omega*TAU)+.5*FUk*sinci(.5*kappa*TAU)*sinci(.5*kappa*TAU)*TAU+.5*FUo*sinc(.5*omega*TAU)*sinc(.5*omega*TAU)*TAU)*TAU-dt;
		if(f0*f1>0) tau0=TAU; else tau1=TAU;
		if(fabs(tau1-tau0)<1e-6*TAU) break;
	}
	if(unlikely(tar)){
		FFF=FF*F;
		u0=(I+F*TAU+1./2.*FF*TAU*TAU+1./6.*FFF*TAU*TAU*TAU)*u0;						//update 4-velocity
		p+=(I+1./2.*F*TAU+1./6.*FF*TAU*TAU+1./24.*FFF*TAU*TAU*TAU)*u0*TAU;				//update 4-position
	}
	else{
		u0=U_k*cosh(kappa*TAU)+F*U_k*sinci(kappa*TAU)*TAU+U_o*cos(omega*TAU)+F*U_o*sinc(omega*TAU)*TAU;	//update 4-velocity
		p+=(U_k*sinci(kappa*TAU)+.5*F*U_k*sinci(.5*kappa*TAU)*sinci(.5*kappa*TAU)*TAU)*TAU+(U_o*sinc(omega*TAU)+.5*F*U_o*sinc(.5*omega*TAU)*sinc(.5*omega*TAU)*TAU)*TAU;	//update 4-position
	}
}
double phase(VecN &u0, VecN &u1){
	double a2=u0.x[0]*u0.x[0]+u0.x[1]*u0.x[1],b2=u1.x[0]*u1.x[0]+u1.x[1]*u1.x[1],c2=(u0.x[0]-u1.x[0])*(u0.x[0]-u1.x[0])+(u0.x[1]-u1.x[1])*(u0.x[1]-u1.x[1]);
	if(u1.x[0]<0) return 2.*M_PI-acos(.5*(a2+b2-c2)/sqrt(a2*b2));
	else return acos(.5*(a2+b2-c2)/sqrt(a2*b2));
}
int main(int argc, char *argv[]){
	double pha,pha0=0.;
	for(dt=1e-8;dt<1.1e-3;dt+=1e-6){
		u0.x[0]=0e3;u0.x[1]=1e3;u0.x[2]=0e3;u0.x[3]=sqrt(1+u0.x[0]*u0.x[0]+u0.x[1]*u0.x[1]+u0.x[2]*u0.x[2]);		//initialize particle parameters
		VecN U0=u0;
		U.x[0]=u0.x[0];U.x[1]=u0.x[1];U.x[2]=u0.x[2];gm=u0.x[3];
#if(defined(a_txt)||defined(c_txt))
		anal_solu();				//analytical solution
#endif
#if(defined(b_txt)||defined(d_txt))
	//	Comp_tau();Update_imp1();		//our algorithm(implementation 1)
		Comp_tau();Update_imp2();		//our algorithm(implementation 2)
	//	Asm_opt();				//optimized
#endif
		pha=phase(U0,u0);
		if(pha0>pha) printf("\n");
		printf("%e\t%e\n",(double)dt,(double)pha);
		fprintf(stderr,"%e\t%e\n",(double)dt,(double)pha);
		pha0=pha;
	}
}
