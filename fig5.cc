#include <stdio.h>
#include <math.h>
#include "mydef.h"
#include "phy_const.h"
#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
//==============constants===============
const int N=4;
#ifdef a_txt
const double dt=1.0e-7;
#endif
#if(defined(b_txt)||defined(d_txt))
const double dt=1.0e-5;
#endif
const double COE0=q_e/m_e/c_0;
//=============declaration==============
struct Mat3;
struct MatN;
//==============structs=================
struct Vec3{						//3-D vector
	double x[3];
	Vec3() {for(int i=0;i<3;i++) x[i]=0;}
	Vec3 operator+(const Vec3 &vec) const {Vec3 temp;for(int i=0;i<3;i++) temp.x[i]=x[i]+vec.x[i];return temp;}
	Vec3 operator-(const Vec3 &vec) const {Vec3 temp;for(int i=0;i<3;i++) temp.x[i]=x[i]-vec.x[i];return temp;}
	double operator*(const Vec3 &vec) const {double temp=0;for(int i=0;i<3;i++) temp+=x[i]*vec.x[i];return temp;}			//dot product
	Vec3 operator*(const double &a) const {Vec3 temp;for(int i=0;i<3;i++) temp.x[i]=x[i]*a;return temp;}				//mult
	friend Vec3 operator*(const double a,const Vec3 &vec) {Vec3 temp;for(int i=0;i<3;i++) temp.x[i]=vec.x[i]*a;return temp;}	//mult
	Mat3 operator&(const Vec3 &vec) const ;		// matrix mult vector
	Mat3 CPM() const ;				// Cross product matrix
};
struct VecN{						//4-D vector
	double x[N];
	VecN() {for(int i=0;i<N;i++) x[i]=0;}
	VecN operator+(const VecN &vec) const {VecN temp;for(int i=0;i<N;i++) temp.x[i]=x[i]+vec.x[i];return temp;}
	VecN operator+=(const VecN &vec) {for(int i=0;i<N;i++) x[i]+=vec.x[i];return *this;}
	VecN operator*(const double &a) const {VecN temp;for(int i=0;i<N;i++) temp.x[i]=x[i]*a;return temp;}				//mult
	friend VecN operator*(const double a,const VecN &vec) {VecN temp;for(int i=0;i<N;i++) temp.x[i]=vec.x[i]*a;return temp;}	//mult
	MatN operator&(const VecN &vec) const ;		// matrix mult vector, defined below;
};
struct Mat3{						//3x3 matrix
	double x[3][3];
	Mat3() {for(int i=0;i<3;i++) for(int j=0;j<3;j++) x[i][j]=0;}
	Mat3 operator*(const double &a) const{Mat3 temp;for(int i=0;i<3;i++) for(int j=0;j<3;j++) temp.x[i][j]=x[i][j]*a;return temp;}	//mult
	Vec3 operator*(const Vec3 &vec) const{Vec3 temp;for(int i=0;i<3;i++) for(int j=0;j<3;j++) temp.x[i]+=x[i][j]*vec.x[j];return temp;}
	friend Vec3 operator*(const Vec3 &a,const Mat3 &mat) {Vec3 temp;for(int i=0;i<3;i++) for(int j=0;j<3;j++) temp.x[i]+=a.x[j]*mat.x[j][i];return temp;}
};
struct MatN{						//4x4 matrix
	double x[N][N];
	MatN() {for(int i=0;i<N;i++) for(int j=0;j<N;j++) x[i][j]=0;}
	MatN(const double x_) {for(int i=0;i<N;i++) {for(int j=0;j<N;j++) x[i][j]=0;x[i][i]=x_;}}
	MatN operator+(const MatN &mat) const {MatN temp;for(int i=0;i<N;i++) for(int j=0;j<N;j++) temp.x[i][j]=x[i][j]+mat.x[i][j];return temp;}
	MatN operator-(const MatN &mat) const {MatN temp;for(int i=0;i<N;i++) for(int j=0;j<N;j++) temp.x[i][j]=x[i][j]-mat.x[i][j];return temp;}
	MatN operator*(const double &a) const{MatN temp;for(int i=0;i<N;i++) for(int j=0;j<N;j++) temp.x[i][j]=x[i][j]*a;return temp;}	//mult
	VecN operator*(const VecN &vec) const{VecN temp;for(int i=0;i<N;i++) for(int j=0;j<N;j++) temp.x[i]+=x[i][j]*vec.x[j];return temp;}
	MatN operator*(const MatN &mat) const{MatN temp;for(int i=0;i<N;i++) for(int j=0;j<N;j++) for(int k=0;k<N;k++) temp.x[i][j]+=x[i][k]*mat.x[k][j];return temp;}
	friend MatN operator*(const double &a,const MatN &mat) {MatN temp;for(int i=0;i<N;i++) for(int j=0;j<N;j++) temp.x[i][j]=mat.x[i][j]*a;return temp;}
	MatN operator/(const double &a) const{MatN temp;for(int i=0;i<N;i++) for(int j=0;j<N;j++) temp.x[i][j]=x[i][j]/a;return temp;}	//div
};
MatN VecN::operator&(const VecN &vec) const{MatN temp;for(int i=0;i<N;i++){temp.x[i][0]=x[i]*vec.x[0];temp.x[i][1]=x[i]*vec.x[1];temp.x[i][2]=x[i]*vec.x[2];temp.x[i][3]=-x[i]*vec.x[3];}return temp;}//matrix mult
Mat3 Vec3::operator&(const Vec3 &vec) const{Mat3 temp;for(int i=0;i<3;i++) for(int j=0;j<3;j++) temp.x[i][j]=x[i]*vec.x[j];return temp;}//matrix mult
Mat3 Vec3::CPM() const{Mat3 mat;mat.x[0][1]=-x[2];mat.x[0][2]=x[1];mat.x[1][0]=x[2];mat.x[1][2]=-x[0];mat.x[2][0]=-x[1];mat.x[2][1]=x[0];return mat;}//Cross product matrix
//===========static variables===========
static double t,tau_,J[3][4],res[3],gm,tau_min,tau_max;	//time, matrix and vector used to solve linear system, and some intermediate variables.
static double coef0,coef1,coef2,coef3,coef4;		//intermediate coefficients
static Vec3 E,B,U;					//electric field, magnetic field, 3-momentum.
static Mat3 IB_inv;
static VecN u0,u1,p,L_U;				//4-velocity in time step n, n+1, 4-position and temporary variable used in iteration.
static MatN F,FF,FFF,L_M,L_M_t;
//============functions=================
inline void FG(){
	F.x[0][0]=      0;F.x[0][1]= B.x[2];F.x[0][2]=-B.x[1];F.x[0][3]=E.x[0];
	F.x[1][0]=-B.x[2];F.x[1][1]=      0;F.x[1][2]= B.x[0];F.x[1][3]=E.x[1];
	F.x[2][0]= B.x[1];F.x[2][1]=-B.x[0];F.x[2][2]=      0;F.x[2][3]=E.x[2];
	F.x[3][0]= E.x[0];F.x[3][1]= E.x[1];F.x[3][2]= E.x[2];F.x[3][3]=     0;
}
inline void Lorentz_matrix(VecN u){			//Lorentz matrix and its inverse, u[4]={gm*beta,gm}
	if(unlikely(u.x[3]==1))
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++) {L_M.x[i][j]=0;L_M_t.x[i][j]=0;}
			L_M.x[i][i]=1;L_M_t.x[i][i]=1;
		}
	else{
		double u2=u.x[3]*u.x[3]-1;		//u2=u*u=gm^2-1
		L_M.x[0][0]=u.x[0]*u.x[0]/u2*(u.x[3]-1)+1;L_M.x[0][1]=u.x[1]*u.x[0]/u2*(u.x[3]-1)  ;L_M.x[0][2]=u.x[2]*u.x[0]/u2*(u.x[3]-1)  ;L_M.x[0][3]=u.x[0];
		L_M.x[1][0]=u.x[0]*u.x[1]/u2*(u.x[3]-1)  ;L_M.x[1][1]=u.x[1]*u.x[1]/u2*(u.x[3]-1)+1;L_M.x[1][2]=u.x[2]*u.x[1]/u2*(u.x[3]-1)  ;L_M.x[1][3]=u.x[1];
		L_M.x[2][0]=u.x[0]*u.x[2]/u2*(u.x[3]-1)  ;L_M.x[2][1]=u.x[1]*u.x[2]/u2*(u.x[3]-1)  ;L_M.x[2][2]=u.x[2]*u.x[2]/u2*(u.x[3]-1)+1;L_M.x[2][3]=u.x[2];
		L_M.x[3][0]=u.x[0]                       ;L_M.x[3][1]=u.x[1]                       ;L_M.x[3][2]=u.x[2]                       ;L_M.x[3][3]=u.x[3];
		L_M_t.x[0][0]= L_M.x[0][0];L_M_t.x[0][1]= L_M.x[0][1];L_M_t.x[0][2]= L_M.x[0][2];L_M_t.x[0][3]=-L_M.x[0][3];
		L_M_t.x[1][0]= L_M.x[1][0];L_M_t.x[1][1]= L_M.x[1][1];L_M_t.x[1][2]= L_M.x[1][2];L_M_t.x[1][3]=-L_M.x[1][3];
		L_M_t.x[2][0]= L_M.x[2][0];L_M_t.x[2][1]= L_M.x[2][1];L_M_t.x[2][2]= L_M.x[2][2];L_M_t.x[2][3]=-L_M.x[2][3];
		L_M_t.x[3][0]=-L_M.x[3][0];L_M_t.x[3][1]=-L_M.x[3][1];L_M_t.x[3][2]=-L_M.x[3][2];L_M_t.x[3][3]= L_M.x[3][3];
	}
}
inline void field(Vec3 &E,Vec3 &B,const double t){	//Electromagnetic fields.
	E.x[0]=0; E.x[1]=0; E.x[2]=0; B.x[0]=0; B.x[1]=0; B.x[2]=-3*COE0;
#ifdef d_txt
	E.x[2]=1*COE0;
#endif
	FG();						//calculate electromagnetic tensor F
	F=L_M_t*F*L_M;					//transform F to computational frame
	E.x[0]=F.x[0][3]; E.x[1]=F.x[1][3]; E.x[2]=F.x[2][3]; B.x[0]=F.x[1][2]; B.x[1]=-F.x[0][2]; B.x[2]=F.x[0][1];	//field in computational frame
}
inline void output(){
	printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",(double)(t),(double)(u0.x[3]),(double)(u0.x[0]/u0.x[3]),(double)(u0.x[1]/u0.x[3]),(double)(u0.x[2]/u0.x[3]),(double)p.x[0],(double)p.x[1],(double)p.x[2]);
	fprintf(stderr,"t:%.4e\tgm:%.4e\tp0:%.4e\tp1:%.4e\tp2:%.4e\n",(double)(t),(double)(u0.x[3]),(double)p.x[0],(double)p.x[1],(double)p.x[2]);
}
inline void Comp_coe(){					//compute coefficients
	coef0=2.*(B*E)*(B*U)+(B*E)*(B*E)*dt,coef1=2.*(U*B.CPM()*E+B*B*gm),coef2=2.*E*U+E*E*dt-B*B*dt,coef3=2.*gm,coef4=-dt;
}
inline double Q_f(double tau){				//quartic function
	return (((coef0*tau+coef1)*tau+coef2)*tau+coef3)*tau+coef4;
}
inline double D_qf(double tau){				//derivative of the quartic function
	return ((4.*coef0*tau+3.*coef1)*tau+2.*coef2)*tau+coef3;
}
inline void Tau_Max(){										//calculate the upper limit of \tau
	double BB=B*B,EE=E*E,BE=B*E,TEMP0=(BB-EE)*(BB-EE),TEMP1=BE*BE,Res;
	if(TEMP0<TEMP1*1e13) Res=sqrt((BB-EE+sqrt((BB-EE)*(BB-EE)+4*BE*BE))/(2*BE*BE));		//when BE*BE is nonzero
	else if(EE>BB) Res=sqrt(1./(EE-BB));							//When BE*BE is approximately zero and EE>BB
	else Res=dt/(1+gm);
	if(Res<dt/(1+gm)) tau_max=Res;else tau_max=dt/(1+gm);					//When EE<=BB, find out tau_max(the upper limit of \tau)
}
inline void Ens_conv(){					//ensure that f is a convex function in [\tau_{min},\tau_{max}]
	if(fabs(coef0)>1e-100){
		double temp=9.*coef1*coef1-24.*coef0*coef2;
		if(temp>=0){
			temp=sqrt(temp);
			double t1=(-3.0*coef1+temp)/(12.*coef0),t2=(-3.0*coef1-temp)/(12.*coef0);
			if(tau_min<t1<tau_max){if(Q_f(t1)>0) tau_max=t1; else tau_min=t1;}
			if(tau_min<t2<tau_max){if(Q_f(t2)>0) tau_max=t2; else tau_min=t2;}
		}
	}
	else if((fabs(coef1)>1e-100)){
		double t1=-coef2/(3.*coef1);
		if(tau_min<t1<tau_max){if(Q_f(t1)>0) tau_max=t1; else tau_min=t1;}
	}
}
inline void Tau_Newt(){		//Newton iteration to solve \tau
	double a=D_qf(tau_min),b=D_qf(tau_max),dt;
	if(a>b){dt=Q_f(tau_min)/a;tau_=tau_min-dt;} else{dt=Q_f(tau_max)/b;tau_=tau_max-dt;}
	while(likely(fabs(dt)>1e-12*tau_)){dt=Q_f(tau_)/D_qf(tau_);tau_-=dt;}
}
inline void Comp_J(){		//compute augmented Jacobian matrix
	J[0][0]=1;		J[0][1]=-B.x[2]*tau_;	J[0][2]=B.x[1]*tau_;	J[0][3]=U.x[0]-B.x[1]*U.x[2]*tau_+B.x[2]*U.x[1]*tau_+E.x[0]*dt;
	J[1][0]=B.x[2]*tau_;	J[1][1]=1;		J[1][2]=-B.x[0]*tau_;	J[1][3]=U.x[1]-B.x[2]*U.x[0]*tau_+B.x[0]*U.x[2]*tau_+E.x[1]*dt;
	J[2][0]=-B.x[1]*tau_;	J[2][1]=B.x[0]*tau_;	J[2][2]=1;		J[2][3]=U.x[2]-B.x[0]*U.x[1]*tau_+B.x[1]*U.x[0]*tau_+E.x[2]*dt;
}
inline void Gauss_Jordan(){	//Gauss Jordan elimination with complete pivoting
	const int N=3;
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
void ZZ(){			//our algorithm(implementation 1)
	Comp_coe();	
	Tau_Max();
	Ens_conv();
	Tau_Newt();
	gm=dt/tau_-gm;
	Comp_J();
	Gauss_Jordan();
	u1=u0;
	u0.x[0]=U.x[0]=res[0];u0.x[1]=U.x[1]=res[1];u0.x[2]=U.x[2]=res[2];gm=u0.x[3]=sqrt(1+u0.x[0]*u0.x[0]+u0.x[1]*u0.x[1]+u0.x[2]*u0.x[2]);
	p+=(u0+u1)*c_0*tau_;
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
inline void Comp_ib_inv(){
	double t2=tau_*tau_;
	IB_inv.x[0][0]=1.+B.x[0]*B.x[0]*t2;		IB_inv.x[0][1]=B.x[2]*tau_+B.x[0]*B.x[1]*t2;	IB_inv.x[0][2]= -B.x[1]*tau_+B.x[0]*B.x[2]*t2;
	IB_inv.x[1][0]= -B.x[2]*tau_+B.x[0]*B.x[1]*t2;	IB_inv.x[1][1]=1+B.x[1]*B.x[1]*t2;		IB_inv.x[1][2]=B.x[0]*tau_+B.x[1]*B.x[2]*t2;
	IB_inv.x[2][0]=B.x[1]*tau_+B.x[0]*B.x[2]*t2;	IB_inv.x[2][1]= -B.x[0]*tau_+B.x[1]*B.x[2]*t2;	IB_inv.x[2][2]=1.+B.x[2]*B.x[2]*t2;
	IB_inv=IB_inv*(1./(1.+B*B*t2));
}
void Comp_tau_opt(){
	field(E,B,t+0.5*dt);				//compute field
	Comp_coe();
	Tau_Max();
	Ens_conv();
	Tau_Newt();
}
void Update_opt(){
	u1=u0;
	Comp_ib_inv();
	U=IB_inv*(2.*U+E*dt)-U;
	u0.x[0]=U.x[0];u0.x[1]=U.x[1];u0.x[2]=U.x[2];u0.x[3]=gm=sqrt(1+U*U);
	p+=(u0+u1)*c_0*tau_;				//update 4-position
}
void ZZ_opt(){			//our algorithm(implementation 2)
	Comp_tau_opt();
	Update_opt();
}
void anal_solu(){
	double TAU,I1=E*E-B*B,I2=E*B;
	double kappa=1/sqrt(2)*sqrt(I1+sqrt(I1*I1+4*I2*I2)),omega=1/sqrt(2)*sqrt(-I1+sqrt(I1*I1+4*I2*I2));
	int tar;
	if(unlikely(kappa<eps_th&&omega<eps_th)){kappa=eps_th;omega=eps_th;tar=1;} else tar=0;
	FF=F*F;
	const MatN I=1.,P_k=(omega*omega*I+FF)/(kappa*kappa+omega*omega),P_o=(kappa*kappa*I-FF)/(kappa*kappa+omega*omega);
	const VecN U_k=P_k*u0,U_o=P_o*u0;
	double tau0=0.,tau1=dt,f0,FUk=F.x[3][0]*U_k.x[0]+F.x[3][1]*U_k.x[1]+F.x[3][2]*U_k.x[2],FUo=F.x[3][0]*U_o.x[0]+F.x[3][1]*U_o.x[1]+F.x[3][2]*U_o.x[2];
	f0=(U_k.x[3]*sinci(kappa*tau0)+U_o.x[3]*sinc(omega*tau0)+.5*FUk*sinci(.5*kappa*tau0)*sinci(.5*kappa*tau0)*tau0+.5*FUo*sinc(.5*omega*tau0)*sinc(.5*omega*tau0)*tau0)*tau0-dt;
/*	for(;;){					//secant method
		double f1=(U_k.x[3]*sinci(kappa*tau1)+U_o.x[3]*sinc(omega*tau1)+.5*FUk*sinci(.5*kappa*tau1)*sinci(.5*kappa*tau1)*tau1+.5*FUo*sinc(.5*omega*tau1)*sinc(.5*omega*tau1)*tau1)*tau1-dt;
		TAU=tau0-f0/(f1-f0)*(tau1-tau0);
		if(fabs(tau1-tau0)<1e-6*TAU) break;
		f0=f1;tau0=tau1;tau1=TAU;
	}*/
	for(;;){					//bisection method
		TAU=.5*(tau1+tau0);
		double f1=(U_k.x[3]*sinci(kappa*TAU)+U_o.x[3]*sinc(omega*TAU)+.5*FUk*sinci(.5*kappa*TAU)*sinci(.5*kappa*TAU)*TAU+.5*FUo*sinc(.5*omega*TAU)*sinc(.5*omega*TAU)*TAU)*TAU-dt;
		if(f0*f1>0) tau0=TAU; else tau1=TAU;
		if(fabs(tau1-tau0)<1e-6*TAU) break;
	}
	if(unlikely(tar)){
		FFF=FF*F;
		u0=(I+F*TAU+1./2.*FF*TAU*TAU+1./6.*FFF*TAU*TAU*TAU)*u0;			//update 4-velocity
		p+=(I+1./2.*F*TAU+1./6.*FF*TAU*TAU+1./24.*FFF*TAU*TAU*TAU)*u0*TAU;	//update 4-position
	}
	else{
		u0=U_k*cosh(kappa*TAU)+F*U_k*sinci(kappa*TAU)*TAU+U_o*cos(omega*TAU)+F*U_o*sinc(omega*TAU)*TAU;	//update 4-velocity
		p+=c_0*((U_k*sinci(kappa*TAU)+.5*F*U_k*sinci(.5*kappa*TAU)*sinci(.5*kappa*TAU)*TAU)*TAU+(U_o*sinc(omega*TAU)+.5*F*U_o*sinc(.5*omega*TAU)*sinc(.5*omega*TAU)*TAU)*TAU);	//update 4-position
	}
}
int main(int argc, char *argv[]){
	L_U.x[1]=0;L_U.x[2]=0;L_U.x[3]=2;L_U.x[0]=sqrt(L_U.x[3]*L_U.x[3]-1);
	Lorentz_matrix(L_U);
	u0.x[0]=0;u0.x[1]=1e3;u0.x[2]=0;u0.x[3]=sqrt(1+u0.x[0]*u0.x[0]+u0.x[1]*u0.x[1]+u0.x[2]*u0.x[2]);		//initialize particle parameters
	u0=L_M_t*u0;					//transform 4-velocity to computational frame
#ifdef a_txt
	for(t=0;t<2.5e-4;t+=dt){
#endif
#if(defined(b_txt)||defined(d_txt))
	for(t=0;t<2e-4;t+=dt){
#endif
		U.x[0]=u0.x[0];U.x[1]=u0.x[1];U.x[2]=u0.x[2];gm=u0.x[3];
		field(E,B,t+0.5*dt);			//compute fields
		u0=L_M*u0;p=L_M*p;			//transform to observational frame
		output();				//output data at time-step n
		u0=L_M_t*u0;p=L_M_t*p;			//back to computational frame
#ifdef a_txt
		anal_solu();				//analytical solution
#endif
#if(defined(b_txt)||defined(d_txt))
	//	ZZ();					//our algorithm(implementation 1)
		ZZ_opt();				//our algorithm(implementation 2)
#endif
	}
	u0=L_M*u0;p=L_M*p;
	output();					//output data of the last time-step
}
