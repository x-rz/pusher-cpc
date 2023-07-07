#include <stdio.h>
#include <math.h>
#include <time.h>
#include "mydef.h"
#include "phy_const.h"
#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
//==============constants===============
const int N=4;
const double dt=1.0e-6;
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
	VecN operator*(const double a) const {VecN temp;for(int i=0;i<N;i++) temp.x[i]=x[i]*a;return temp;}				//mult
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
volatile double vlt=1.0,vlt1;				//variables used to ensure the compiler does not ignore some code due to optimization
static double t,tau_,J[4][5],res[4],gm,tau_min,tau_max;	//time, matrix and vector used to solve linear system, and some intermediate variables
static double coef0,coef1,coef2,coef3,coef4;		//intermediate coefficients
static Vec3 E,B,U;					//electric field, magnetic field, 3-momentum
static Mat3 IB_inv;
static VecN u0,u1,p;					//4-velocity in time step n, n+1, 4-position
static MatN F,FF,FFF;
//============functions=================
inline void FG(){
	F.x[0][0]=      0;F.x[0][1]= B.x[2];F.x[0][2]=-B.x[1];F.x[0][3]=E.x[0];
	F.x[1][0]=-B.x[2];F.x[1][1]=      0;F.x[1][2]= B.x[0];F.x[1][3]=E.x[1];
	F.x[2][0]= B.x[1];F.x[2][1]=-B.x[0];F.x[2][2]=      0;F.x[2][3]=E.x[2];
	F.x[3][0]= E.x[0];F.x[3][1]= E.x[1];F.x[3][2]= E.x[2];F.x[3][3]=     0;
}
inline void field(Vec3 &E,Vec3 &B,const double t){	//Electromagnetic fields. vlt used to ensure the compiler does not ignore related code
	E.x[0]=0*COE0*vlt; E.x[1]=5.196152*COE0*vlt; E.x[2]=2*COE0*vlt; B.x[0]=0*COE0*vlt; B.x[1]=1.732051*COE0*vlt; B.x[2]=-6*COE0*vlt;
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
	J[0][0]=1;J[0][1]=-B.x[2]*tau_;J[0][2]=B.x[1]*tau_;J[0][3]=-E.x[0]*tau_;J[0][4]=u0.x[0]-B.x[1]*u0.x[2]*tau_+B.x[2]*u0.x[1]*tau_+E.x[0]*u0.x[3]*tau_;
	J[1][0]=B.x[2]*tau_;J[1][1]=1;J[1][2]=-B.x[0]*tau_;J[1][3]=-E.x[1]*tau_;J[1][4]=u0.x[1]-B.x[2]*u0.x[0]*tau_+B.x[0]*u0.x[2]*tau_+E.x[1]*u0.x[3]*tau_;
	J[2][0]=-B.x[1]*tau_;J[2][1]=B.x[0]*tau_;J[2][2]=1;J[2][3]=-E.x[2]*tau_;J[2][4]=u0.x[2]-B.x[0]*u0.x[1]*tau_+B.x[1]*u0.x[0]*tau_+E.x[2]*u0.x[3]*tau_;
	J[3][0]=-E.x[0]*tau_;J[3][1]=-E.x[1]*tau_;J[3][2]=-E.x[2]*tau_;J[3][3]=1;J[3][4]=u0.x[3]+E.x[0]*u0.x[0]*tau_+E.x[1]*u0.x[1]*tau_+E.x[2]*u0.x[2]*tau_;
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
void Update(){
	gm=dt/tau_-gm;
	Comp_J();
	Gauss_Jordan();
	u1=u0;
	u0.x[0]=U.x[0]=res[0];u0.x[1]=U.x[1]=res[1];u0.x[2]=U.x[2]=res[2];gm=u0.x[3]=res[3];
	p+=(u0+u1)*c_0*tau_;				//update 4-position
}
inline void Comp_ib_inv(){
	double t2=tau_*tau_;
	IB_inv.x[0][0]=1.+B.x[0]*B.x[0]*t2;		IB_inv.x[0][1]=B.x[2]*tau_+B.x[0]*B.x[1]*t2;	IB_inv.x[0][2]= -B.x[1]*tau_+B.x[0]*B.x[2]*t2;
	IB_inv.x[1][0]= -B.x[2]*tau_+B.x[0]*B.x[1]*t2;	IB_inv.x[1][1]=1+B.x[1]*B.x[1]*t2;		IB_inv.x[1][2]=B.x[0]*tau_+B.x[1]*B.x[2]*t2;
	IB_inv.x[2][0]=B.x[1]*tau_+B.x[0]*B.x[2]*t2;	IB_inv.x[2][1]= -B.x[0]*tau_+B.x[1]*B.x[2]*t2;	IB_inv.x[2][2]=1.+B.x[2]*B.x[2]*t2;
	IB_inv=IB_inv*(1./(1.+B*B*t2));
}
inline void Tau_Newt_opt(){	//Newton iteration to solve \tau
	double temp=.5*(tau_min+tau_max);			//bisection
	if(Q_f(temp)>0) tau_max=temp;else tau_min=temp;		//bisection
	temp=.5*(tau_min+tau_max);				//bisection
	if(Q_f(temp)>0) tau_max=temp;else tau_min=temp;		//bisection
	double a=D_qf(tau_min),b=D_qf(tau_max),dt;
	if(a>b){dt=Q_f(tau_min)/a;tau_=tau_min-dt;} else{dt=Q_f(tau_max)/b;tau_=tau_max-dt;}
	while(likely(fabs(dt)>1e-12*tau_)){dt=Q_f(tau_)/D_qf(tau_);tau_-=dt;}
}
void Comp_tau_opt(){
	field(E,B,t+0.5*dt);				//compute field
	Comp_coe();
	Tau_Max();
	Ens_conv();
	Tau_Newt_opt();
}
void Update_opt(){
	u1=u0;
	Comp_ib_inv();
	U=IB_inv*(2.*U+E*dt)-U;
	u0.x[0]=U.x[0];u0.x[1]=U.x[1];u0.x[2]=U.x[2];u0.x[3]=gm=sqrt(1+U*U);
	p+=(u0+u1)*c_0*tau_;				//update 4-position
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
double TAU,I1,I2,kappa,omega,tau0,tau1,f0,f1,FUk,FUo;
int tar;
const MatN I=1.;
MatN P_k,P_o;
VecN U_k,U_o;
void Comp_tau_as(){
	field(E,B,t+0.5*dt);				//compute field
	I1=E*E-B*B,I2=E*B;
	kappa=1/sqrt(2)*sqrt(I1+sqrt(I1*I1+4*I2*I2)),omega=1/sqrt(2)*sqrt(-I1+sqrt(I1*I1+4*I2*I2));
	if(unlikely(kappa<eps_th&&omega<eps_th)){kappa=eps_th;omega=eps_th;tar=1;}
	else tar=0;
	FG();FF=F*F;
	P_k=(omega*omega*I+FF)/(kappa*kappa+omega*omega),P_o=(kappa*kappa*I-FF)/(kappa*kappa+omega*omega);
	U_k=P_k*u0,U_o=P_o*u0;
	tau0=0.,tau1=dt; FUk=F.x[3][0]*U_k.x[0]+F.x[3][1]*U_k.x[1]+F.x[3][2]*U_k.x[2],FUo=F.x[3][0]*U_o.x[0]+F.x[3][1]*U_o.x[1]+F.x[3][2]*U_o.x[2];
	f0=(U_k.x[3]*sinci(kappa*tau0)+U_o.x[3]*sinc(omega*tau0)+.5*FUk*sinci(.5*kappa*tau0)*sinci(.5*kappa*tau0)*tau0+.5*FUo*sinc(.5*omega*tau0)*sinc(.5*omega*tau0)*tau0)*tau0-dt;
	for(;;){					//secant method	
		f1=(U_k.x[3]*sinci(kappa*tau1)+U_o.x[3]*sinc(omega*tau1)+.5*FUk*sinci(.5*kappa*tau1)*sinci(.5*kappa*tau1)*tau1+.5*FUo*sinc(.5*omega*tau1)*sinc(.5*omega*tau1)*tau1)*tau1-dt;
		TAU=tau0-f0/(f1-f0)*(tau1-tau0);
		if(unlikely(tau1-tau0<1e-12*tau1)) break;
		f0=f1;tau0=tau1;tau1=TAU;
	}
}
void Update_as(){
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
void ZZ(){
	Comp_tau();
	Update();
}
void ZZ_opt(){
	Comp_tau_opt();
	Update_opt();
}
void anal_solu(){
	Comp_tau_as();
	Update_as();
}
void sample_error(double *sa, int sn, double &aver, double &sd){
	aver=0;	for(int i=0;i<sn;i++) aver+=sa[i]; aver/=sn;			//average
	double temp=0; for(int i=0;i<sn;i++) temp+=((sa[i]-aver)*(sa[i]-aver));	//errors square sum
	sd=sqrt(temp/(sn-1));							//sample standard deviation 
} 
int main(int argc, char *argv[]){
	int sn=5,stt,end;			//sample number, start clock, end clock
	double dta[sn],aver,sd;			//data buffer, average, SD
	double tau[4],coef; int n;
	u0.x[0]=-1.732052e+03;u0.x[1]=0.000000e+00;u0.x[2]=1.000000e+03;u0.x[3]=2.000001e+03;		//initialize particle parameters
	U.x[0]=u0.x[0];U.x[1]=u0.x[1];U.x[2]=u0.x[2];gm=u0.x[3];
	for(int i=0;i<sn;i++){
		stt=clock();
		for(t=0;t<=1e7*dt;t+=dt)
#ifdef anal_t
			Comp_tau_as();vlt1=tau_;		//vlt1 used to ensure the compiler does not ignore related code
#elif(defined(anal_p))
			Update_as();vlt1=p.x[0];		//vlt1 used to ensure the compiler does not ignore related code
#elif(defined(imp1_t))
			Comp_tau();vlt1=tau_;			//vlt1 used to ensure the compiler does not ignore related code
#elif(defined(imp1_p))
			Update();vlt1=p.x[0];			//vlt1 used to ensure the compiler does not ignore related code
#elif(defined(imp2_t))
			Comp_tau_opt();vlt1=tau_;		//vlt1 used to ensure the compiler does not ignore related code
#elif(defined(imp2_p))
			Update_opt();vlt1=u0.x[0];vlt1=p.x[0];	//vlt1 used to ensure the compiler does not ignore related code
#endif		
		end=clock();
		dta[i]=1e9*(end-stt)/CLOCKS_PER_SEC/1e7;	//in unit of nano second
	}
	sample_error(dta, sn, aver, sd);
	printf("%.4e\t%.4e",aver,sd);
}
