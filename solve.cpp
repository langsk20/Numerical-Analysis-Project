#include <bits/stdc++.h>
using namespace std;
const int K=10;
const int N=1e6;
const int M=1e4;
const double pi=acos(-1);
const double eps=1e-10;
const double C=4.0;
const double TOL=1e-6; 
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
struct treenode {
	double alphal,alphar,betal,betar,deltal,deltar,mul,mur;
	double l,r;
	bool is_leaf;
	int lc,rc,fa;
}T[N*4+10];
int tot=1;
double node[K+1];
void Simple_Solve(double l,double r,int x) {
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
	s1=Mul(Cf,s1);
	s2=Mul(Cf,s2);
	double sum1=0,sum2=0;
	for (int i=0;i<K;i+=2) {
		sum1+=2.0/(1-i*i)*s1.v[i+1];
		sum2+=2.0/(1-i*i)*s2.v[i+1];
	}
	T[x].deltal=sum1;
	T[x].deltar=sum2;
	for (int i=1;i<=K;i++) sigma.v[i]=psi_l(node[i]);
	sigma=Linear_Solver(P,sigma);
	for (int i=1;i<=K;i++) {
		s1.v[i]=sigma.v[i]*gl(node[i])*(r-l)/2;
		s2.v[i]=sigma.v[i]*gr(node[i])*(r-l)/2; 
	}
	s1=Mul(Cf,s1);
	s2=Mul(Cf,s2);
	sum1=sum2=0;
	for (int i=0;i<K;i+=2) {
		sum1+=2.0/(1-i*i)*s1.v[i+1];
		sum2+=2.0/(1-i*i)*s2.v[i+1];
	}
	T[x].alphal=sum1;
	T[x].alphar=sum2;
	for (int i=1;i<=K;i++) sigma.v[i]=psi_r(node[i]);
	sigma=Linear_Solver(P,sigma);
	for (int i=1;i<=K;i++) {
		s1.v[i]=sigma.v[i]*gl(node[i])*(r-l)/2;
		s2.v[i]=sigma.v[i]*gr(node[i])*(r-l)/2; 
	}
	s1=Mul(Cf,s1);
	s2=Mul(Cf,s2);
	sum1=sum2=0;
	for (int i=0;i<K;i+=2) {
		sum1+=2.0/(1-i*i)*s1.v[i+1];
		sum2+=2.0/(1-i*i)*s2.v[i+1];
	}
	T[x].betal=sum1;
	T[x].betar=sum2;
	return;
}
void Update(int x,int lc,int rc) {
	double Delta=1-T[rc].alphar*T[lc].betal;
	T[x].alphal=(1-T[rc].alphal)*(T[lc].alphal-T[lc].betal*T[rc].alphar)/Delta+T[rc].alphal;
	T[x].alphar=T[rc].alphar*(1-T[lc].betar)*(1-T[lc].alphal)/Delta+T[lc].alphar;
	T[x].betal=T[lc].betal*(1-T[rc].betar)*(1-T[rc].alphal)/Delta+T[rc].betal;
	T[x].betar=(1-T[lc].betar)*(T[rc].betar-T[lc].betal*T[rc].alphar)/Delta+T[lc].betar;
	T[x].deltal=(1-T[rc].alphal)/Delta*T[lc].deltal+T[rc].deltal+
		(T[rc].alphal-1)*T[lc].betal/Delta*T[rc].deltar;
	T[x].deltar=(1-T[lc].betar)/Delta*T[rc].deltar+T[lc].deltar+
		(T[lc].betar-1)*T[rc].alphar/Delta*T[lc].deltal;
}
void Down(int x,int lc,int rc) {
	T[lc].mul=T[x].mul;
	T[rc].mur=T[x].mur;
	matrix S;
	S.n=2;
	S.a[1][1]=S.a[2][2]=1;
	S.a[1][2]=T[rc].alphar;
	S.a[2][1]=T[lc].betal;
	vec t;
	t.n=2;
	t.v[1]=T[x].mur*(1-T[rc].betar)-T[rc].deltar;
	t.v[2]=T[x].mul*(1-T[lc].alphal)-T[lc].deltal;
	t=Linear_Solver(S,t);
	T[lc].mur=t.v[1];
	T[rc].mul=t.v[2];
}
vector<int> qu,qu2;
void Downward(int x) {
	if (T[x].is_leaf) {
		Simple_Solve(T[x].l,T[x].r,x);
		return;
	}
	Downward(T[x].lc);
	Downward(T[x].rc);
	Update(x,T[x].lc,T[x].rc);
}
void Upward(int x) {
	if (x==1) return;
	Upward(T[x].fa);
	Down(T[x].fa,T[T[x].fa].lc,T[T[x].fa].rc); 
}
double Jl[N+10],Jr[N+10];
vec poly[N+10];
double Si[N+10];
void Calc(int i0,int x,double jl,double jr) {
	double l=T[x].l,r=T[x].r;
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
	vec sigma1,sigma2;
	sigma1.n=sigma2.n=K;
	for (int i=1;i<=K;i++) sigma1.v[i]=psi_l(node[i])*T[x].mul;
	for (int i=1;i<=K;i++) sigma2.v[i]=psi_r(node[i])*T[x].mur;
	sigma1=Linear_Solver(P,sigma1);
	sigma2=Linear_Solver(P,sigma2);
	for (int i=1;i<=K;i++) sigma.v[i]+=sigma1.v[i]+sigma2.v[i];
	vec si=Mul(Cf,sigma);
	Si[i0]=abs(si.v[K-1])+abs(si.v[K-2]+si.v[K]);
	vec s1=sigma,s2=sigma;
	for (int i=1;i<=K;i++) {
		s1.v[i]*=gl(node[i])*(r-l)/2;
		s2.v[i]*=gr(node[i])*(r-l)/2; 
	}
	s1=Mul(Il,s1);
	s2=Mul(Ir,s2);
	poly[i0].n=K;
	for (int i=1;i<=K;i++) {
		double x=node[i];
		double u=Ui(x);
		u+=gr(x)/s(x)*(jl+s1.v[i]);
		u+=gl(x)/s(x)*(jr+s2.v[i]);
		poly[i0].v[i]=u;
	}
}
void Get() {
	int sz=qu.size();
	for (int i=0;i<sz;i++) T[qu[i]].is_leaf=1;
	Downward(1);
	for (int i=0;i<sz;i++) T[qu[i]].is_leaf=0;
	for (int i=0;i<sz;i++) Upward(qu[i]);
	Jl[0]=0;
	for (int i=0;i<sz;i++)
		Jl[i+1]=Jl[i]+T[qu[i]].deltal+T[qu[i]].mul*T[qu[i]].alphal+T[qu[i]].mur*T[qu[i]].betal;
	Jr[sz-1]=0;
	for (int i=sz-1;i>=1;i--) 
		Jr[i-1]=Jr[i]+T[qu[i]].deltar+T[qu[i]].mul*T[qu[i]].alphar+T[qu[i]].mur*T[qu[i]].betar;
	for (int i=0;i<sz;i++)
		Calc(i,qu[i],Jl[i],Jr[i]);
}
double result[M+1],nxt[M+1];
void Output() {
	int sz=qu.size();
	int k=0;
	for (int i=0;i<sz;i++) poly[i]=Mul(Cf,poly[i]);
	for (int i=0;i<=M;i++) {
		double now=L+(R-L)/M*i;
		if (T[qu[k]].r+eps<now) k++;
		double l=T[qu[k]].l,r=T[qu[k]].r;
		result[i]=0;
		for (int j=1;j<=K;j++)
			result[i]+=poly[k].v[j]*cos((j-1)*acos(2*(now-l)/(r-l)-1)); 
	}
}
void Build(int x) {
	double mid=(T[x].l+T[x].r)/2;
	T[x].lc=++tot;
	T[tot].fa=x;
	T[tot].l=T[x].l;
	T[tot].r=mid;
	T[x].rc=++tot;
	T[tot].fa=x;
	T[tot].l=mid;
	T[tot].r=T[x].r; 
}
void Solve() {
	T[1].fa=0;
	T[1].mul=T[1].mur=0;
	T[1].l=L;
	T[1].r=R;
	qu.clear();
	qu.push_back(1);
	bool judge=false;
	while (1) {
		cerr<<"one loop"<<endl;
		Get();
		Output();
		double sum1=0,sum2=0;
		for (int i=0;i<=M;i++) {
			sum1+=(result[i]-nxt[i])*(result[i]-nxt[i]);
			sum2+=(result[i]+nxt[i])*(result[i]+nxt[i]);
		}
		double test=sqrt(sum1/sum2);
		for (int i=0;i<=M;i++) nxt[i]=result[i];
		if (test>TOL) {
			double Sdiv=0;
			for (int i=0;i<(int)qu.size();i++) Sdiv=max(Sdiv,Si[i]);
			Sdiv/=pow(2,C);
			qu2.clear();
			for (int i=0;i<(int)qu.size();i++) {
				if (Si[i]>=Sdiv) {
					Build(qu[i]);
					qu2.push_back(T[qu[i]].lc);
					qu2.push_back(T[qu[i]].rc);
				}
				else {
					if (i<(int)qu.size()-1&&T[qu[i]].fa==T[qu[i+1]].fa
						&&Si[i]+Si[i+1]<Sdiv/pow(2,K)) {
							qu2.push_back(T[qu[i]].fa);
							i++; 
					}
					else qu2.push_back(qu[i]);
				}
			}
			qu=qu2;
		}
		else {
			if (!judge) {
				qu2.clear();
				for (int i=0;i<(int)qu.size();i++) {
					Build(qu[i]);
					qu2.push_back(T[qu[i]].lc);
					qu2.push_back(T[qu[i]].rc);
				}
				qu=qu2;
				judge=true;
			}
			else break;
		}
	}
	for (int i=0;i<=M;i++) {
		printf("%lf %lf\n",L+(R-L)/M*i,result[i]);
	}
}
int main() {
	freopen("out.txt","w",stdout);
	L=-pi;
	R=pi;
	zeta_l0=1;zeta_l1=1;
	zeta_r0=2;zeta_r1=3;
	Gamma_l=0;Gamma_r=1;
	ans.clear(); 
	Pre();
	Solve();
	return 0;
}
