#include <immintrin.h>
#ifndef LIB_H
#define LIB_H
union Mat3;
union MatN;
union __attribute__((aligned(32))) Vec3{						//3-D vector
	double x[4];
	__m256d v;
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) Vec3(){for(int i=0;i<4;i++) x[i]=0;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) Vec3 operator+(const Vec3 &vec)const{Vec3 temp;for(int i=0;i<3;i++) temp.x[i]=x[i]+vec.x[i];return temp;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) Vec3 operator-(const Vec3 &vec)const{Vec3 temp;for(int i=0;i<3;i++) temp.x[i]=x[i]-vec.x[i];return temp;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) double operator*(const Vec3 &vec)const{double temp=0;for(int i=0;i<3;i++) temp+=x[i]*vec.x[i];return temp;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) Vec3 operator*(const double &a)const{Vec3 temp;for(int i=0;i<3;i++) temp.x[i]=x[i]*a;return temp;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) friend Vec3 operator*(const double &a,const Vec3 &vec){Vec3 temp;for(int i=0;i<3;i++) temp.x[i]=vec.x[i]*a;return temp;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) Mat3 CPM() const; // Cross product matrix
};
union __attribute__((aligned(32))) VecN{						//4-D vector
	double x[4];
	__m256d v;
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) VecN(){for(int i=0;i<4;i++) x[i]=0;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) VecN operator+(const VecN &vec)const{VecN temp;for(int i=0;i<4;i++) temp.x[i]=x[i]+vec.x[i];return temp;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) VecN operator+=(const VecN &vec){for(int i=0;i<4;i++) x[i]+=vec.x[i];return *this;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) VecN operator*(const double &a)const{VecN temp;for(int i=0;i<4;i++) temp.x[i]=x[i]*a;return temp;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) friend VecN operator*(const double &a,const VecN &vec){VecN temp;for(int i=0;i<4;i++) temp.x[i]=vec.x[i]*a;return temp;}
};
union __attribute__((aligned(32))) Mat3{						//3x3 matrix
	double x[4][4];
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) Mat3(){for(int i=0;i<3;i++) for(int j=0;j<3;j++) x[i][j]=0;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) Mat3 operator*(const double &a)const{Mat3 temp;for(int i=0;i<3;i++) for(int j=0;j<3;j++) temp.x[i][j]=x[i][j]*a;return temp;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) Vec3 operator*(const Vec3 &vec)const{Vec3 temp;for(int i=0;i<3;i++) for(int j=0;j<3;j++) temp.x[i]+=x[i][j]*vec.x[j];return temp;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) friend Vec3 operator*(const Vec3 &vec,const Mat3 &mat){Vec3 temp;for(int i=0;i<3;i++) for(int j=0;j<3;j++) temp.x[i]+=vec.x[j]*mat.x[j][i];return temp;}
};
union __attribute__((aligned(32))) MatN{						//4x4 matrix
	double x[4][4];
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) MatN(){for(int i=0;i<4;i++) for(int j=0;j<4;j++) x[i][j]=0;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) MatN(const double x_){for(int i=0;i<4;i++) {for(int j=0;j<4;j++) x[i][j]=0;x[i][i]=x_;}}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) MatN operator+(const MatN &mat)const{MatN temp;for(int i=0;i<4;i++) for(int j=0;j<4;j++) temp.x[i][j]=x[i][j]+mat.x[i][j];return temp;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) MatN operator-(const MatN &mat)const{MatN temp;for(int i=0;i<4;i++) for(int j=0;j<4;j++) temp.x[i][j]=x[i][j]-mat.x[i][j];return temp;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) MatN operator*(const double &a)const{MatN temp;for(int i=0;i<4;i++) for(int j=0;j<4;j++) temp.x[i][j]=x[i][j]*a;return temp;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) VecN operator*(const VecN &vec)const{VecN temp;for(int i=0;i<4;i++) for(int j=0;j<4;j++) temp.x[i]+=x[i][j]*vec.x[j];return temp;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) MatN operator*(const MatN &mat)const{MatN temp;for(int i=0;i<4;i++) for(int j=0;j<4;j++) for(int k=0;k<4;k++) temp.x[i][j]+=x[i][k]*mat.x[k][j];return temp;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) friend MatN operator*(const double &a,const MatN &mat){MatN temp;for(int i=0;i<4;i++) for(int j=0;j<4;j++) temp.x[i][j]=mat.x[i][j]*a;return temp;}
	__inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) MatN operator/(const double &a)const{MatN temp;for(int i=0;i<4;i++) for(int j=0;j<4;j++) temp.x[i][j]=x[i][j]/a;return temp;}
};

Mat3 Vec3::CPM() const {Mat3 mat;mat.x[0][1]=-x[2];mat.x[0][2]=x[1];mat.x[1][0]=x[2];mat.x[1][2]=-x[0];mat.x[2][0]=-x[1];mat.x[2][1]=x[0];return mat;}
#endif
