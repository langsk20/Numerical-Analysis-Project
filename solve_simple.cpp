#include <bits/stdc++.h>
using namespace std;
const int K=50;
const double pi=acos(-1);
const double eps=1e-10;
const double C=4.0;
double L,R,zeta_l0,zeta_l1,zeta_r0,zeta_r1,Gamma_l,Gamma_r;
int q0;
struct UI {
	double k,b;
}ui;
struct vec {
	int n;
	double v[K+1];
	void print() {
		printf("%d\n",n);
		for (int i=1;i<=n;i++) printf("%lf ",v[i]);
		printf("\n");
	}
};
struct matrix {
	int n;
	double a[K+1][K+1];
	void print() {
		printf("%d\n",n);
		for (int i=1;i<=n;i++) {
			for (int j=1;j<=n;j++) printf("%lf ",a[i][j]);
			printf("\n");
		}
	}
}Cf,Cb,Ils,Irs,Il,Ir;
vector <pair<double,double>> ans;
void Addpoint(double x,double u) {
	ans.push_back(make_pair(x,u));
}
double gl(double x) {
	return (q0==0?(zeta_l0*(x-L)-zeta_l1):(zeta_l1*cosh(x-L)-zeta_l0*sinh(x-L)));
}
double gr(double x) {
	return (q0==0?(zeta_r0*(x-R)-zeta_r1):(zeta_r1*cosh(x-R)-zeta_r0*sinh(x-R)));
}
double D_gl(double x) {
	return (q0==0?zeta_l0:(zeta_l1*sinh(x-L)-zeta_l0*cosh(x-L)));
}
double D_gr(double x) {
	return (q0==0?zeta_r0:(zeta_r1*sinh(x-R)-zeta_r0*cosh(x-R)));
}
double p(double x) {
	return 0;
}
double q(double x) {
	return 1;
}
double f(double x) {
	return 0;
}
double p_tilde(double x) {
	return p(x);
}
double q_tilde(double x) {
	return q(x)-q0;
}
double s(double x) {
	return gl(x)*D_gr(x)-D_gl(x)*gr(x);
}
double psi_l(double x) {
	return (p_tilde(x)*D_gr(x)+q_tilde(x)*gr(x))/s(x);
}
double psi_r(double x) {
	return (p_tilde(x)*D_gl(x)+q_tilde(x)*gl(x))/s(x);
}
double Ui(double x) {
	return ui.k*x+ui.b;
}
double D_Ui(double x) {
	return ui.k;
}
double D2_Ui(double x) {
	return 0;
} 
double f_tilde(double x) {
	return f(x)-D2_Ui(x)-p(x)*D_Ui(x)-q(x)*Ui(x);
}
matrix Multi(matrix A,matrix B) {
	matrix c;
	c.n=A.n;
	int n=c.n;
	for (int i=1;i<=n;i++)
		for (int j=1;j<=n;j++) {
			c.a[i][j]=0;
			for (int k=1;k<=n;k++)
				c.a[i][j]+=A.a[i][k]*B.a[k][j]; 
		}
	return c;
}
vec Mul(matrix A,vec x) {
	vec b;
	b.n=A.n;
	int n=b.n;
	for (int i=1;i<=n;i++) {
		b.v[i]=0;
		for (int j=1;j<=n;j++)
			b.v[i]+=A.a[i][j]*x.v[j];
	}
	return b;
}
vec Linear_Solver(matrix A,vec b) { //Gauss Elimination
	int n=A.n;
	for (int i=1;i<=n;i++) {
		double max0=0;
		int line;
		for (int j=i;j<=n;j++) 
		  if (abs(A.a[j][i])>max0) {
			  max0=abs(A.a[j][i]);
			  line=j;
		  }
		if (abs(max0)<=eps) {
			printf("Wrong\n");
			exit(0);
		}
		if (line!=i) {
		  for (int j=1;j<=n;j++) swap(A.a[i][j],A.a[line][j]);
		  swap(b.v[i],b.v[line]);
		}
		for (int j=1;j<=n;j++) {
			if (j==i) continue;
			double t=A.a[j][i]/A.a[i][i];
			for (int k=i;k<=n;k++) A.a[j][k]-=t*A.a[i][k];
			b.v[j]-=t*b.v[i];
		}
	}
	vec ans;
	ans.n=n;
	for (int i=1;i<=n;i++) 
		ans.v[i]=b.v[i]/A.a[i][i];
	return ans;
}
void Pre() {
	if (abs(zeta_l0)>=abs(zeta_l1)||abs(zeta_r0)>=abs(zeta_r1)) q0=0;
	else q0=-1;
	ui.k=(Gamma_l*zeta_r0-Gamma_r*zeta_l0)/
		(zeta_l0*zeta_r0*(L-R)+zeta_r0*zeta_l1-zeta_l0*zeta_r1);
	if (abs(zeta_l0)>eps) ui.b=(Gamma_l-zeta_l1*ui.k-zeta_l0*L*ui.k)/zeta_l0;
	else ui.b=(Gamma_r-zeta_r1*ui.k-zeta_r0*R*ui.k)/zeta_r0;
	Cf.n=K;
	for (int i=1;i<=Cf.n;i++) Cf.a[1][i]=1.0/K;
	for (int i=2;i<=Cf.n;i++)
		for (int j=1;j<=Cf.n;j++)
			Cf.a[i][j]=2.0/K*cos((i-1)*(2*K-2*j+1)*pi/2/K);
	Cb.n=K;
	for (int i=1;i<=Cb.n;i++)
		for (int j=1;j<=Cb.n;j++)
			Cb.a[i][j]=cos((j-1)*(2*K-2*i+1)*pi/2/K);
	Ils.n=K;
	Ils.a[2][1]=1;
	Ils.a[2][3]=-0.5;
	for (int i=2;i<K;i++) {
		Ils.a[i+1][i]=0.5/i;
		if (i+1<K) Ils.a[i+1][i+2]=-0.5/i;
	}
	for (int i=2,sgn=1;i<=K;i++,sgn=-sgn)
		for (int j=1;j<=K;j++) Ils.a[1][j]+=sgn*Ils.a[i][j];
	Irs.n=K;
	Irs.a[2][1]=-1;
	Irs.a[2][3]=0.5;
	for (int i=2;i<K;i++) {
		Irs.a[i+1][i]=-0.5/i;
		if (i+1<K) Irs.a[i+1][i+2]=0.5/i;
	}
	for (int i=2;i<=K;i++)
		for (int j=1;j<=K;j++) Irs.a[1][j]-=Irs.a[i][j];
	Il=Multi(Cb,Multi(Ils,Cf));
	Ir=Multi(Cb,Multi(Irs,Cf));
}
double node[K+1];
void Simple_Solve(double l,double r) {
	matrix P;
	P.n=K;
	for (int i=1;i<=K;i++) node[i]=(r-l)/2*cos((2*K-2*i+1)*pi/2/K)+(l+r)/2;
	for (int i=1;i<=K;i++)
		for (int j=1;j<=K;j++) {
			double now=0;
			if (i==j) now=1;
			now+=psi_l(node[i])*gl(node[j])*Il.a[i][j]*(r-l)/2;
			now+=psi_r(node[i])*gr(node[j])*Ir.a[i][j]*(r-l)/2;
			P.a[i][j]=now;
		}
	vec sigma;
	sigma.n=K;
	for (int i=1;i<=K;i++) sigma.v[i]=f_tilde(node[i]);
	sigma=Linear_Solver(P,sigma);
	vec s1=sigma,s2=sigma;
	for (int i=1;i<=K;i++) {
		s1.v[i]*=gl(node[i])*(r-l)/2;
		s2.v[i]*=gr(node[i])*(r-l)/2; 
	}
	s1=Mul(Il,s1);
	s2=Mul(Ir,s2);
	for (int i=1;i<=K;i++) {
		double x=node[i];
		double u=Ui(x);
		u+=gr(x)/s(x)*s1.v[i];
		u+=gl(x)/s(x)*s2.v[i];
		Addpoint(x,u);
	}
}
void Print() {
	for (int i=0;i<(int)ans.size();i++)
		printf("%lf %lf\n",ans[i].first,ans[i].second);
}
int main() {
	freopen("output.txt","w",stdout);
	L=-pi;
	R=pi;
	zeta_l0=1;zeta_l1=1;
	zeta_r0=2;zeta_r1=3;
	Gamma_l=0;Gamma_r=1;
	ans.clear(); 
	Pre();
	Simple_Solve(L,R);
	Print();
	return 0;
}
