#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include <omp.h>

#define N 100

void cargar_int(int * A, char * ruta, unsigned int X_size, unsigned int Y_size);
void cargar_double(double * A, char * ruta, unsigned int X_size, unsigned int Y_size);
unsigned int idx(unsigned int i, unsigned int j, unsigned int stride);
int min_gradPr(double * x, int * Gamma, int G_col,int * Delta,double * c_a,int * reg_A);
void proy(double * h, int * Gamma, int G_col);
double T(double * x,double * df,int * Delta,double * c_a,int * reg_A,int G_col,int flag);
int reg_armijoPr(double * x,double f,double * df,int G_col,int * Gamma,int * Delta,double * c_a,int * reg_A);
void qp_i(double * x0, double * h_aux, int m);
void qp_e(double * p_k, double * g_k, int * W_k, int m, double * lambda);

unsigned int c_a_x=9006,c_a_y=1,D_x=355561,D_y=2,G_x=2990,G_y=1;
int maxit=50,maxit_r=50,nn=0;
double tol,s=100;

int main(int argc, char **argv) {
	assert(2==argc);
	D_x = atof(argv[1]);
	size_t c_a_size = c_a_x * c_a_y * sizeof(double);
	size_t D_size = D_x * D_y * sizeof(int);
	size_t G_size = G_x * G_y * sizeof(int);
	int reg_A[4],salida=1;
	double * c_a = malloc(c_a_size);
	int * Delta = malloc(D_size);
	int * Gamma = malloc(G_size);
	cargar_int(reg_A,"min/data/input/reg_A",1,4); //reg_A=load(’data/input/reg_A’);
	cargar_double(c_a,"min/data/input/costo_aristas",c_a_x,c_a_y); //c_a=load(’data/input/costo_aristas’);
	cargar_int(Delta,"grafo_completo/data/output/Delta",D_x,D_y); //Delta=load(’data/input/Delta’);
	cargar_int(Gamma,"grafo_completo/data/output/Gamma",G_x,G_y); //Gamma=load(’data/input/Gamma’);
	int G_col=0;
	for (int i=0;i<G_x;i++){
		G_col+=Gamma[i];
	}
	size_t h_size = G_col * sizeof(double);
	tol=DBL_EPSILON*G_col;
	double * h = malloc(h_size);
	int acum=0,k=0;
	for (int i=0;i<G_col;i++){	//for i=1:rows(Gamma)
		if (i==acum){			//h0(min(find(Gamma(i,:))))=1;
			h[i]=1;				//endfor
			acum+=Gamma[k];		//if (length(h0)<columns(Gamma))
			k++;				//h0(columns(Gamma))=0;
		} else {				//endif
			h[i]=0;				//h0=h0(:);
		}
	}
//[h,f,err,k,salida]=min_gradPr(@T,h0,@(x)qp_s(Gamma,x),[100 1e-12],@reg_armijoPr,[100 1000 1 1 0.5 1e -4]);
	int flag=0;
	double start=omp_get_wtime();
	while(salida && nn<500){
		if (flag){
			s=s/2;
			if (s<DBL_EPSILON*10){
				s=1000;
				flag=0;
				salida=1;
				continue;
			}
			salida=min_gradPr(h,Gamma,G_col,Delta,c_a,reg_A);
		} else if (salida==3){
			s=s/2;
			salida=min_gradPr(h,Gamma,G_col,Delta,c_a,reg_A);
			flag=1;
		} else {
			s=s*2;
			salida=min_gradPr(h,Gamma,G_col,Delta,c_a,reg_A);
		}
	}
	printf("%f\n",omp_get_wtime()-start);
	double * v = malloc(c_a_x * sizeof(double));
	for (int i=0;i<c_a_x;i++){
		v[i]=0;
	}
	for (int i=0;i<D_x;i++){
		v[-1+Delta[idx(i,1,D_y)]]+=h[-1+Delta[idx(i,0,D_y)]]; // v=Delta'*h;
	}
	FILE * v_file;
	v_file = fopen ("min/data/output/v", "w+");
	for (int i=0;i<c_a_x;i++){
		fprintf(v_file,"%lf\n",v[i]);
	}
	fclose(v_file);
	free(c_a);
	free(Delta);
	free(Gamma);
	free(h);
	free(v);
	return 0;
}

void cargar_int(int * A, char * ruta, unsigned int X_size, unsigned int Y_size){
	FILE * archivo;
	archivo = fopen (ruta, "r");
	for (int i=0;i<X_size;i++){
		for (int j=0;j<Y_size;j++){
			fscanf(archivo,"%d",&A[idx(i,j,Y_size)]);
		}
	}
	fclose(archivo);
}

void cargar_double(double * A, char * ruta, unsigned int X_size, unsigned int Y_size){
	FILE * archivo;
	archivo = fopen (ruta, "r");
	for (int i=0;i<X_size;i++){
		for (int j=0;j<Y_size;j++){
			fscanf(archivo,"%lf",&A[idx(i,j,Y_size)]);
		}
	}
	fclose(archivo);
}

unsigned int idx(unsigned int x, unsigned int y, unsigned int stride) {
    return x * stride + y;
}

int min_gradPr(double * x, int * Gamma, int G_col,int * Delta,double * c_a,int * reg_A){
	size_t h_size = G_col * sizeof(double);
	double f, err; //,fold;
	double * df = malloc(h_size);
	double * h = malloc(h_size);
	int ki;
	proy(x,Gamma,G_col); // x=PrX(x);
	f=T(x,df,Delta,c_a,reg_A,G_col,1); // [f,df]=fun(x,[1,1]);
	for (int i=0;i<G_col;i++){
		h[i]=x[i]-df[i];
	}
	proy(h,Gamma,G_col);
	err=0;
	for (int i=0; i<G_col;i++){
		if (fabs(x[i]-h[i])> (10*tol*D_x*D_x)/(G_col*c_a_x*G_col) ) {
			err += (x[i] - h[i])*(x[i] - h[i]); // err=norm(x-PrX(x-df));
		}
	}
	err=sqrt(err);
	for (int i=0;i<maxit;i++){ // for k=0:maxit
		printf("iter_min=%d fun=%.16lf err=%.16lf s=%lf\n",nn++,f,err,s); // printf ("iter_min=%d fun=%e err=%e\n",k,f,err)
		if (err<tol){ // if err<tol
			free(df);
			free(h);
			return 0; // salida=0;
		} // end
		ki=reg_armijoPr(x,f,df,G_col,Gamma,Delta,c_a,reg_A); // [x,t,ki]=regla(fun,x,PrX,f,df,pr);
		if (ki==-1){
			free(df);
			free(h);
			return 4; // salida=3;
		}
		if (ki>=maxit_r){ // if ki>=pr(1)
			free(df);
			free(h);
			return 3; // salida=3;
		} // end
		// fold=f;
		f=T(x,df,Delta,c_a,reg_A,G_col,1); // [f,df]=fun(x,[1,1]);
		/*
		if (fold==f){
			free(df);
			free(h);
			return 2; // salida=0;
		}
		*/
		for (int i=0;i<G_col;i++){
			h[i]=x[i]-df[i];
		}
		proy(h,Gamma,G_col);
		err=0;
		for (int i=0; i<G_col;i++){
			if (fabs(x[i]-h[i])> (10*tol*D_x*D_x)/(G_col*c_a_x*G_col) ){
				err += (x[i] - h[i])*(x[i] - h[i]); // err=norm(x-PrX(x-df));
			}
		}
		err=sqrt(err);
	} // end
	free(df);
	free(h);
	return 1;
}

/*
	min_gradPr(@T,h0,@(x)qp_s(Gamma,x),[100 1e-12],@reg_armijoPr,[100 1000 1 1 0.5 1e -4]);
function [x,f,err,k,salida]=min_gradPr(fun,x,PrX,parmet,funregla,parregla)
% ejecutar [x,f,err,k,salida]=min_cauchy(fun,x0,PrX)
%        o [x,f,err,k,salida]=min_cauchy(fun,x0,PrX,parmet)
%        o [x,f,err,k,salida]=min_cauchy(fun,x0,PrX,parmte,funregla)
%        o [x,f,err,k,salida]=min_cauchy(fun,x0,PrX,parmte,funregla,parregla)
%
% donde parmet(1)=cantidad maxima de iteraciones (default=100)
%       parmet(2)=tolerancia de error del gradiente (default=1e-7)
%       funregla=funcion regla de busqueda lineal (default=@reg_armijo)
%       parregla=parametros para la regla escogida
%
% salida=0, gradiente menor que tolerancia.
% salida=1, maximo de iteraciones alcanzado.
% salida=3, maximo de iteraciones en busqueda lineal.
pm=[100,1e-7]; regla=@reg_armijoPr; pr=100;
if nargin>3, pm(1:length(parmet))=parmet; end
if nargin>4&&~isempty(funregla), regla=funregla; end
if nargin==6, pr(1:length(parregla))=parregla; end
maxit=pm(1); tol=pm(2);
*/

double T(double * h,double * df,int * Delta,double * c_a,int * reg_A,int G_col,int flag){
	size_t v_size = c_a_x * sizeof(double);
	double * v = malloc(v_size);
	double * w = malloc(v_size);
	double f;
	for (int i=0;i<c_a_x;i++){
		v[i]=0;
	}
	for (int i=0;i<D_x;i++){
		v[-1+Delta[idx(i,1,D_y)]]+=h[-1+Delta[idx(i,0,D_y)]]; // v=Delta'*h;
	}
	memcpy(w,v,v_size); // w=v;
	for (int i=0;i<reg_A[0];i++){
		w[i]=v[i]-exp(-v[i]); // w(1:reg_A(1))=v(1:reg_A(1))-exp(-v(1:reg_A(1)));
	}
	f=0;
	for (int i=0;i<c_a_x;i++){
		f+=c_a[i]*w[i]; // f=sum(c_a.*w);
	}
	if (flag){
		for (int i=0;i<c_a_x;i++){
			w[i]=1;
		}
		for (int i=0;i<reg_A[0];i++){
			w[i]=1+exp(-v[i]); // z=[1+exp(-v(1:reg_A(1)));ones(length(v)-reg_A(1),1)];
		}
		for (int i=0;i<c_a_x;i++){
			w[i]=c_a[i]*w[i]; // z=c_a.*z;
		}
		for (int i=0;i<G_col;i++){
			df[i]=0;
		}
		for (int i=0;i<D_x;i++){
			df[-1+Delta[idx(i,0,D_y)]]+=w[-1+Delta[idx(i,1,D_y)]]; // Df=Delta*z;
		}
	}
	free(v);
	free(w);
	return f;
}

int reg_armijoPr(double * x,double f,double * df,int G_col,int * Gamma,int * Delta,double * c_a,int * reg_A){
	size_t h_size = G_col * sizeof(double);
	double * xs = malloc(h_size);
	double * d = malloc(h_size);
	double * xn = malloc(h_size);
	double alpha=1,fn,d_aux=0,s_aux=s;
	int k=0;
	for (int i=0;i<G_col;i++){
		xs[i]=x[i]-s_aux*df[i];
	}
	proy(xs,Gamma,G_col); // xs=PrX(x-s*df); 
	for (int i=0;i<G_col;i++){
		d[i]=xs[i]-x[i]; // d=xs-x;
	}
	for (int i=0;i<G_col;i++){
		xn[i]=x[i]+alpha*d[i]; // xn=x+alpha*d;
	}
	fn=T(xn,df,Delta,c_a,reg_A,G_col,0); // fn=fun(xn);
	for (int i=0;i<G_col;i++){
		d_aux+=df[i]*d[i];
	}
	if (d_aux>=0){
		free(xs);
		free(d);
		free(xn);
		return -1;
	}
//	printf("d_aux=%.16f ",d_aux);
	while (fn>f+alpha*(1e-4)*d_aux && k<maxit_r){ // while fn>f+alpha*sig*df'*d && k<=maxit
		//printf("    iter_reg=%d fun_n=%.14lf fun=%.14lf\n",k+1,fn,f); //printf ("iter_reg=%d fun_n=%e fun=%e\n",k,fn,f)
		s_aux=0.5*s_aux;
//		alpha=0.5*alpha; // alpha=ra*alpha;
		for (int i=0;i<G_col;i++){
			xs[i]=x[i]-s_aux*df[i];
		}
		proy(xs,Gamma,G_col); // xs=PrX(x-s*df); 
		for (int i=0;i<G_col;i++){
			d[i]=xs[i]-x[i]; // d=xs-x;
		}
		for (int i=0;i<G_col;i++){
			xn[i]=x[i]+alpha*d[i]; // xn=x+alpha*d;
		}
		fn=T(xn,df,Delta,c_a,reg_A,G_col,0); // fn=fun(xn);
		k++; // k=k+1;
	} // end
	if (k<maxit_r){
		for (int i=0;i<G_col;i++){
			x[i]=xn[i];
		}
	}	
//	proy(x,Gamma,G_col);
	free(xs);
	free(d);
	free(xn);
	return k;
}

void proy(double * h, int * Gamma, int G_col){
	size_t i_size = G_x * sizeof(int);
	double h_aux[N],x0[N];
	int * i_aux = malloc(i_size);
	int j=0,m; // j=1;
	for (int i =0;i<G_x;i++){
		i_aux[i]=j;
		j+=Gamma[i];
	}
	#pragma omp parallel shared(Gamma,i_aux,h) private(m,h_aux,x0) //num_threads(n_hilos)
	{
	#pragma omp for schedule(guided)
	for (int i=0;i<G_x;i++){ // for i=1:rows(Gamma)
		m=Gamma[i]; // m=length(find(Gamma(i,:)));
		for (int k=i_aux[i];k<i_aux[i]+m;k++){
			h_aux[k-i_aux[i]]=-h[k]; // h_aux=h0([j:j+m-1]);
		}
		for (int k=0;k<m;k++){
			x0[k]=1.0/((double) m); // x0=(1/m)*ones(m,1);
		}
		qp_i(x0,h_aux,m);
		for (int k=i_aux[i];k<i_aux[i]+m;k++){
			h[k]=x0[k-i_aux[i]]; // x1([j:j+m-1])=qp(x0,speye(m),-h_aux,ones(1,m),1,zeros(m,1),[]);
		}
	} // endfor
	}
	free(i_aux);
}

/*
min q(x) = 0.5 x^T G x + x^T c

s.a. 	a_i^T x =  b_i 	i in E
		a_i^T x >= b_i	i in I

G=Id
c=-h_aux
a_i=ones for i in E={1}
b_i=1	 for i in E={1}
a_i=e_i  for i in I={1,...,m}
b_i=0    for i in I={1,...,m}
*/

void qp_i(double * x0, double * h_aux, int m){
	// Compute a feasible starting point x_0
	int * W = malloc((m+1)*sizeof(int));
	double g_k[N],p_k[N],lambda[N+1],lambda_min,alpha;
	int flag,j,block;
	// Set W_0 to be a subset of the active constraints at x_0
	for (int i=0;i<m;i++){
		if (x0[i]==0){
			W[i]=1;
		} else {
			W[i]=0;
		}
	}
	W[m]=1;
	for (int k=0;k<100;k++){ // for k = 0, 1, 2, . . .
		// g_k = x_k +c
		for (int i=0;i<m;i++){
			g_k[i]=x0[i]+h_aux[i];
		}
		// p_k = min  0.5 p^T p + g_k^T p
		//		 s.a.       p_i = 0, i in W_k-{m}
		//			  sum_i p_i = 1 
		qp_e(p_k,g_k,W,m,lambda);
		flag=1;
		for (int i=0;i<m && flag;i++){
			if (p_k[i] != 0){
				flag=0;
			}
		}
		if (flag){ // if p_k=0
			// Compute Lagrange multipliers λ_i that satisfy (16.42),
			// sum_{i in W_k} a_i λ_i = g = G x_k + c 
			flag=1;
			for (int i=0;i<m && flag;i++){
				if (W[i]==1 && lambda[i]<0){
					flag=0;
				}
			}
			if (flag){ // if λ_i ≥ 0 for all i ∈ W_k ∩ I
				// stop with solution x* = x_k;
				free(W);
				return;
			} else { // else
				// j <- argmin_{j in W_k ∩ I} λ_i
				lambda_min=1e12;
				j=-1;
				for (int i=0;i<m;i++){
					if (W[i]==1 && lambda[i]<lambda_min){
						lambda_min=lambda[i];
						j=i;
					}
				}
				assert(j!=-1);
				// x_{k+1} <- x_k 
				W[j]=0; // W_k <- W_k \ {j}
			}
		} else { // else
			block=-1;
			// α_k = min{i not in W_k , p_k_i <0} (- x_k_i)/(p_k_i)
			alpha=1e12;
			for (int i=0;i<m;i++){
				if (W[i]==0 && p_k[i]<0 && alpha > -x0[i]/p_k[i]) {
					alpha=-x0[i]/p_k[i];
					block=i;
				}
			}
			// α_k = min (1 , α_k )
			if (alpha > 1){
				alpha=1;
				block=-1;
			}
			// x_{k+1} <- x_k + α_k * p_k
			for (int i=0;i<m;i++){
				x0[i]+= alpha * p_k[i];
			}
			if (block!=-1){ // if there are blocking constraints
				// Obtain W_{k+1} by adding one of the blocking
				// constraints to W_k
				W[block]=1;
			}
			// else
			//     W_{k+1} <- W_k
		}
	} // end (for)
	free(W);
}

/*
p_k = min  0.5 p^T p + g_k^T p
	  s.a.       p_i = 0, i in W_k-{m}
		   sum_i p_i = 0 
 
(    Id    -Id_{W_k}^T -ones ) ( p_k )   ( -g_k )
( Id_{W_k}       0       0  ) (  λ* ) = (   0  )	
(   ones         0       0  ) ( λ_m )   (   0  )

p_1 - lambda_1 - lambda_m+1 = - g_k_1
p_2 - lambda_2 - lambda_m+1 = - g_k_2
...
p_m - lambda_m - lambda_m+1 = - g_k_m
p_1 = 0
p_2 = 0
...
p_m = 0
p_1 + p_2 + ... + p_m = 0
 
i in W_k =>
p_i = 0
lambda_i + lambda_m+1 = g_k_i 

i not in W_k =>
p_i - lambda_m+1 = - g_k_i

p_1 + p_2 + ... + p_m = 1 =>
L * lambda_m+1 - sum_{i not in W_k} g_k_i = 0

=>
 
lambda_m+1 = ( sum_{i not in W_k} g_k_i ) / L

=> tengo todo lo que quiero

*/

void qp_e(double * p_k, double * g_k, int * W_k, int m, double * lambda){
	int L=0;
	double g=0;
	for (int i=0;i<m;i++){
		if (W_k[i]==0){
			L++;
			g+=g_k[i];
		}
	}
	lambda[m]=g/L;
	for (int i=0;i<m;i++){
		if (W_k[i]){
			p_k[i]=0;
			lambda[i]=g_k[i]-lambda[m];
		} else {
			p_k[i]=-g_k[i]+lambda[m];
		}
	}
}
