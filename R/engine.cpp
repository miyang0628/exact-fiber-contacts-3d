// engine.cpp -- Rcpp kernels for 3D fiber contact counts (R port of engine.py)
// Contact: axis distance <= 2r and contact point (midpoint of closest points)
// inside the counting box C=[0,a1]x[0,a2]x[0,a3].
#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <cmath>
using namespace Rcpp;

// closest points of P(s)=p1+s d1, Q(t)=p2+t d2, s,t in [0,1] (Ericson 2005, 5.1.9)
static inline void seg_closest(const double* p1,const double* d1,const double* p2,const double* d2,double& s,double& t){
  double r0=p1[0]-p2[0], r1=p1[1]-p2[1], r2=p1[2]-p2[2];
  double a=d1[0]*d1[0]+d1[1]*d1[1]+d1[2]*d1[2];
  double e=d2[0]*d2[0]+d2[1]*d2[1]+d2[2]*d2[2];
  double f=d2[0]*r0+d2[1]*r1+d2[2]*r2;
  double c=d1[0]*r0+d1[1]*r1+d1[2]*r2;
  double b=d1[0]*d2[0]+d1[1]*d2[1]+d1[2]*d2[2];
  double den=a*e-b*b;
  if(den>1e-14*a*e){ s=(b*f-c*e)/den; s=std::min(1.0,std::max(0.0,s)); } else s=0.0;
  t=(b*s+f)/e;
  if(t<0.0){ t=0.0; s=std::min(1.0,std::max(0.0,-c/a)); }
  else if(t>1.0){ t=1.0; s=std::min(1.0,std::max(0.0,(b-c)/a)); }
}

// fiber i: X[i,], T[i,], L[i]; fiber j shifted by sh (torus image)
static inline int pair_contact(const NumericMatrix& X,const NumericMatrix& T,const NumericVector& L,
                               int i,int j,const double* sh,double r,const double* a,bool gate){
  double p1[3],d1[3],p2[3],d2[3];
  for(int k=0;k<3;k++){
    d1[k]=L[i]*T(i,k); p1[k]=X(i,k)-0.5*d1[k];
    d2[k]=L[j]*T(j,k); p2[k]=X(j,k)+sh[k]-0.5*d2[k];
  }
  double s,t; seg_closest(p1,d1,p2,d2,s,t);
  double dd=0, pc[3];
  for(int k=0;k<3;k++){ double u=p1[k]+s*d1[k], v=p2[k]+t*d2[k]; dd+=(u-v)*(u-v); pc[k]=0.5*(u+v); }
  if(dd>4.0*r*r) return 0;
  if(gate) for(int k=0;k<3;k++) if(pc[k]<0.0||pc[k]>a[k]) return 0;
  return 1;
}

static void aabb(const NumericMatrix& X,const NumericMatrix& T,const NumericVector& L,double r,
                 std::vector<double>& lo,std::vector<double>& hi){
  int n=X.nrow(); lo.resize(3*n); hi.resize(3*n);
  for(int i=0;i<n;i++) for(int k=0;k<3;k++){
    double h=0.5*L[i]*std::fabs(T(i,k)); lo[3*i+k]=X(i,k)-h-r; hi[3*i+k]=X(i,k)+h+r; }
}

//' Sweep-and-prune contact list (0-based indices i<j) in box C
// [[Rcpp::export]]
IntegerMatrix pairs_box(NumericMatrix X,NumericMatrix T,NumericVector L,double r,NumericVector box){
  int n=X.nrow(); std::vector<double> lo,hi; aabb(X,T,L,r,lo,hi);
  std::vector<int> ord(n); std::iota(ord.begin(),ord.end(),0);
  std::sort(ord.begin(),ord.end(),[&](int u,int v){return lo[3*u]<lo[3*v];});
  double a[3]={box[0],box[1],box[2]}, sh[3]={0,0,0};
  std::vector<int> I,J;
  for(int ii=0;ii<n;ii++){ int i=ord[ii];
    for(int jj=ii+1;jj<n;jj++){ int j=ord[jj];
      if(lo[3*j]>hi[3*i]) break;
      if(lo[3*j+1]>hi[3*i+1]||lo[3*i+1]>hi[3*j+1]) continue;
      if(lo[3*j+2]>hi[3*i+2]||lo[3*i+2]>hi[3*j+2]) continue;
      if(pair_contact(X,T,L,i,j,sh,r,a,true)){ I.push_back(std::min(i,j)); J.push_back(std::max(i,j)); }
    }}
  IntegerMatrix P(I.size(),2);
  for(size_t m=0;m<I.size();m++){ P(m,0)=I[m]; P(m,1)=J[m]; }
  return P;
}

//' Nested counts N_n for each n in ngrid (prefix samples), cube/box
// [[Rcpp::export]]
NumericVector count_box(NumericMatrix X,NumericMatrix T,NumericVector L,double r,NumericVector box,IntegerVector ngrid){
  IntegerMatrix P=pairs_box(X,T,L,r,box); int n=X.nrow();
  std::vector<double> hits(n,0.0);
  for(int m=0;m<P.nrow();m++) hits[P(m,1)]+=1.0;
  std::partial_sum(hits.begin(),hits.end(),hits.begin());
  NumericVector out(ngrid.size()); for(int g=0;g<ngrid.size();g++) out[g]=hits[ngrid[g]-1];
  return out;
}

//' Nested counts on the unit torus (minimum image of the midpoint difference, no gate)
// [[Rcpp::export]]
NumericVector count_torus(NumericMatrix X,NumericMatrix T,NumericVector L,double r,IntegerVector ngrid){
  int n=X.nrow(); std::vector<double> hits(n,0.0); double a[3]={1,1,1};
  for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){
    double sh[3]; bool ok=true;
    for(int k=0;k<3;k++){ double dx=X(j,k)-X(i,k); double s=-std::floor(dx+0.5); sh[k]=s;
      if(std::fabs(dx+s)>0.5*(L[i]*std::fabs(T(i,k))+L[j]*std::fabs(T(j,k)))+2*r){ok=false;break;} }
    if(ok && pair_contact(X,T,L,i,j,sh,r,a,false)) hits[j]+=1.0;
  }
  std::partial_sum(hits.begin(),hits.end(),hits.begin());
  NumericVector out(ngrid.size()); for(int g=0;g<ngrid.size();g++) out[g]=hits[ngrid[g]-1];
  return out;
}

//' Brute force count (all pairs)
// [[Rcpp::export]]
double count_brute(NumericMatrix X,NumericMatrix T,NumericVector L,double r,NumericVector box){
  int n=X.nrow(); double a[3]={box[0],box[1],box[2]}, sh[3]={0,0,0}; double c=0;
  for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) c+=pair_contact(X,T,L,i,j,sh,r,a,true);
  return c;
}

//' Sweep-and-prune: contacts and candidate pairs
// [[Rcpp::export]]
NumericVector count_sap(NumericMatrix X,NumericMatrix T,NumericVector L,double r,NumericVector box){
  int n=X.nrow(); std::vector<double> lo,hi; aabb(X,T,L,r,lo,hi);
  std::vector<int> ord(n); std::iota(ord.begin(),ord.end(),0);
  std::sort(ord.begin(),ord.end(),[&](int u,int v){return lo[3*u]<lo[3*v];});
  double a[3]={box[0],box[1],box[2]}, sh[3]={0,0,0}; double c=0,cand=0;
  for(int ii=0;ii<n;ii++){ int i=ord[ii];
    for(int jj=ii+1;jj<n;jj++){ int j=ord[jj];
      if(lo[3*j]>hi[3*i]) break;
      if(lo[3*j+1]>hi[3*i+1]||lo[3*i+1]>hi[3*j+1]) continue;
      if(lo[3*j+2]>hi[3*i+2]||lo[3*i+2]>hi[3*j+2]) continue;
      cand++; c+=pair_contact(X,T,L,i,j,sh,r,a,true); }}
  return NumericVector::create(_["contacts"]=c,_["candidates"]=cand);
}

//' Exact cell list (Algorithm 1): contacts, candidates, registrations
// [[Rcpp::export]]
NumericVector count_cell(NumericMatrix X,NumericMatrix T,NumericVector L,double r,NumericVector box,double h,double org){
  int n=X.nrow(); std::vector<double> lo,hi; aabb(X,T,L,r,lo,hi);
  double mx=*std::max_element(hi.begin(),hi.end());
  long ng=(long)std::floor((mx-org)/h)+1, nc=ng*ng*ng;
  std::vector<long> cnt(nc+1,0); long nreg=0;
  auto cell=[&](double v){return (long)std::floor((v-org)/h);};
  for(int i=0;i<n;i++){
    long g0[3],g1[3]; for(int k=0;k<3;k++){g0[k]=cell(lo[3*i+k]); g1[k]=cell(hi[3*i+k]);}
    for(long x=g0[0];x<=g1[0];x++) for(long y=g0[1];y<=g1[1];y++) for(long z=g0[2];z<=g1[2];z++){ cnt[(x*ng+y)*ng+z+1]++; nreg++; }
  }
  std::partial_sum(cnt.begin(),cnt.end(),cnt.begin());
  std::vector<long> fill(cnt.begin(),cnt.end()-1); std::vector<int> items(nreg);
  for(int i=0;i<n;i++){
    long g0[3],g1[3]; for(int k=0;k<3;k++){g0[k]=cell(lo[3*i+k]); g1[k]=cell(hi[3*i+k]);}
    for(long x=g0[0];x<=g1[0];x++) for(long y=g0[1];y<=g1[1];y++) for(long z=g0[2];z<=g1[2];z++){ long c=(x*ng+y)*ng+z; items[fill[c]++]=i; }
  }
  double a[3]={box[0],box[1],box[2]}, sh[3]={0,0,0}; double contacts=0,cand=0;
  for(long c=0;c<nc;c++){
    long cz=c%ng, cy=(c/ng)%ng, cx=c/(ng*ng); long cc[3]={cx,cy,cz};
    for(long p=cnt[c];p<cnt[c+1];p++){ int i=items[p];
      for(long q=p+1;q<cnt[c+1];q++){ int j=items[q]; bool ok=true, own=true;
        for(int k=0;k<3;k++){ double l=std::max(lo[3*i+k],lo[3*j+k]), u=std::min(hi[3*i+k],hi[3*j+k]);
          if(l>u){ok=false;break;} if(cell(l)!=cc[k]) own=false; }
        if(!ok||!own) continue;
        cand++; contacts+=pair_contact(X,T,L,i,j,sh,r,a,true); }}
  }
  return NumericVector::create(_["contacts"]=contacts,_["candidates"]=cand,_["registrations"]=(double)nreg);
}

//' Unbiased estimate of |S| h_{1,r}(omega1) by sampling partners in a box containing K_r
//' (partner direction: isotropic with z-component scaled by c, then normalized). Uses R's RNG.
// [[Rcpp::export]]
NumericVector h1_est(NumericVector X1,NumericVector T1,double l1,NumericVector lens,NumericVector cum,
                     double r,int M,NumericVector box,double c){
  NumericMatrix Xa(2,3),Ta(2,3); NumericVector La(2);
  for(int k=0;k<3;k++){Xa(0,k)=X1[k]; Ta(0,k)=T1[k];} La[0]=l1;
  double a[3]={box[0],box[1],box[2]}, sh[3]={0,0,0}; double s1=0,s2=0;
  for(int m=0;m<M;m++){
    double v[3]={R::norm_rand(),R::norm_rand(),R::norm_rand()*c};
    double nv=std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); for(int k=0;k<3;k++) v[k]/=nv;
    double u=R::unif_rand(); int idx=0; while(u>cum[idx]) idx++; double l2=lens[idx];
    double cg=T1[0]*v[0]+T1[1]*v[1]+T1[2]*v[2];
    double e2[3]; for(int k=0;k<3;k++) e2[k]=v[k]-cg*T1[k];
    double sg=std::sqrt(e2[0]*e2[0]+e2[1]*e2[1]+e2[2]*e2[2]); if(sg<1e-12) continue;
    for(int k=0;k<3;k++) e2[k]/=sg;
    double e3[3]={T1[1]*e2[2]-T1[2]*e2[1], T1[2]*e2[0]-T1[0]*e2[2], T1[0]*e2[1]-T1[1]*e2[0]};
    double A=0.5*(l1+l2*std::fabs(cg))+2*r, B=0.5*l2*sg+2*r, Cc=2*r, Vb=8*A*B*Cc;
    double x=(2*R::unif_rand()-1)*A, y=(2*R::unif_rand()-1)*B, z=(2*R::unif_rand()-1)*Cc;
    for(int k=0;k<3;k++){ Xa(1,k)=X1[k]+x*T1[k]+y*e2[k]+z*e3[k]; Ta(1,k)=v[k]; } La[1]=l2;
    double val=Vb*pair_contact(Xa,Ta,La,0,1,sh,r,a,true); s1+=val; s2+=val*val;
  }
  double mean=s1/M; return NumericVector::create(mean,(s2/M-mean*mean)/M);
}
