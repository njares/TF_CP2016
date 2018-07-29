#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

void cargar_int(int * A, char * ruta, unsigned int X_size, unsigned int Y_size);
void cargar_float(float * A, char * ruta, unsigned int X_size, unsigned int Y_size);
unsigned int idx(unsigned int i, unsigned int j, unsigned int stride);
int min_gradPr(float * x, int * Gamma, int G_col,int * Delta,float * c_a,int * reg_A);
void qp_s(float * h, int * Gamma, int G_col);
float T(float * x,float * df,int * Delta,float * c_a,int * reg_A,int G_col,int flag);
int reg_armijoPr(float * x,float f,float * df,int G_col,int * Gamma,int * Delta,float * c_a,int * reg_A);

unsigned int c_a_x=9006,c_a_y=1,D_x=355746,D_y=2,G_x=2990,G_y=1;
int maxit=100,maxit_r=100;
float tol=1e-12,s=1;

int main(){
	size_t c_a_size = c_a_x * c_a_y * sizeof(float);
	size_t D_size = D_x * D_y * sizeof(int);
	size_t G_size = G_x * G_y * sizeof(int);
	int reg_A[4],salida;
	float * c_a = malloc(c_a_size);
	int * Delta = malloc(D_size);
	int * Gamma = malloc(G_size);
	cargar_int(reg_A,"data/input/reg_A",1,4); //reg_A=load(’data/input/reg_A’);
	cargar_float(c_a,"data/input/costo_aristas",c_a_x,c_a_y); //c_a=load(’data/input/costo_aristas’);
	cargar_int(Delta,"data/input/Delta",D_x,D_y); //Delta=load(’data/input/Delta’);
	cargar_int(Gamma,"data/input/Gamma",G_x,G_y); //Gamma=load(’data/input/Gamma’);
	int G_col=0;
	for (int i=0;i<G_x;i++){
		G_col+=Gamma[i];
	}
	size_t h_size = G_col * sizeof(float);
	float * h = malloc(h_size);
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
	salida=min_gradPr(h,Gamma,G_col,Delta,c_a,reg_A);
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

void cargar_float(float * A, char * ruta, unsigned int X_size, unsigned int Y_size){
	FILE * archivo;
	archivo = fopen (ruta, "r");
	for (int i=0;i<X_size;i++){
		for (int j=0;j<Y_size;j++){
			fscanf(archivo,"%f",&A[idx(i,j,Y_size)]);
		}
	}
	fclose(archivo);
}

unsigned int idx(unsigned int x, unsigned int y, unsigned int stride) {
    return x * stride + y;
}

int min_gradPr(float * x, int * Gamma, int G_col,int * Delta,float * c_a,int * reg_A){
	size_t h_size = G_col * sizeof(float);
	float f, err;
	float * df = malloc(h_size);
	float * h = malloc(h_size);
	int ki;
	qp_s(x,Gamma,G_col); // x=PrX(x);
	f=T(x,df,Delta,c_a,reg_A,G_col,1); // [f,df]=fun(x,[1,1]);
	for (int i=0;i<G_col;i++){
		h[i]=x[i]-df[i];
	}
	qp_s(h,Gamma,G_col);
	err=0;
	for (int i=0; i<=G_col;i++){
		err += (x[i] - h[i])*(x[i] - h[i]); // err=norm(x-PrX(x-df));
	}
	err=sqrt(err);
	for (int i=0;i<maxit;i++){ // for k=0:maxit
		printf("iter_min=%d fun=%f err=%f\n",i+1,f,err); // printf ("iter_min=%d fun=%e err=%e\n",k,f,err)
		if (err<tol){ // if err<tol
			return 0; // salida=0;
		} // end
		ki=reg_armijoPr(x,f,df,G_col,Gamma,Delta,c_a,reg_A); // [x,t,ki]=regla(fun,x,PrX,f,df,pr);
		if (ki>=maxit_r){ // if ki>=pr(1)
			return 3; // salida=3;
		} // end
		f=T(x,df,Delta,c_a,reg_A,G_col,1); // [f,df]=fun(x,[1,1]);
		for (int i=0;i<G_col;i++){
			h[i]=x[i]-df[i];
		}
		qp_s(h,Gamma,G_col);
		err=0;
		for (int i=0; i<=G_col;i++){
			err += (x[i] - h[i])*(x[i] - h[i]); // err=norm(x-PrX(x-df));
		}
		err=sqrt(err);		
	} // end
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

void qp_s(float * h, int * Gamma, int G_col){
	size_t h_size = G_col * sizeof(float);
	float h_aux[100];
	int j=0,m; // j=1;
	for (int i=0;i<G_x;i++){ // for i=1:rows(Gamma)
		m=Gamma[i]; // m=length(find(Gamma(i,:)));
		for (int k=j;k<j+m;k++){
			h_aux[k-j]=h[k]; // h_aux=h0([j:j+m-1]);
		}
/*
		x0=(1/m)*ones(m,1);
		x1([j:j+m-1])=qp(x0,speye(m),-h_aux,ones(1,m),1,zeros(m,1),[]);
*/
		j=j+m; // j=j+m;	
	} // endfor
}

float T(float * h,float * df,int * Delta,float * c_a,int * reg_A,int G_col,int flag){
	size_t v_size = c_a_x * sizeof(float);
	float * v = malloc(v_size);
	float * w = malloc(v_size);
	float f;
	for (int i=0;i<c_a_x;i++){
		v[i]=0;
	}
	for (int i=0;i<D_x;i++){
		v[Delta[idx(i,1,D_y)]]+=h[Delta[idx(i,0,D_y)]]; // v=Delta'*h;
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
			df[Delta[idx(i,0,D_y)]]+=w[Delta[idx(i,1,D_y)]]; // Df=Delta*z;
		}
	}
	return f;
}

int reg_armijoPr(float * x,float f,float * df,int G_col,int * Gamma,int * Delta,float * c_a,int * reg_A){
	size_t h_size = G_col * sizeof(float);
	float * xs = malloc(h_size);
	float * d = malloc(h_size);
	float * xn = malloc(h_size);
	float alpha=1,fn,d_aux=0;
	int k=0;
	for (int i=0;i<G_col;i++){
		xs[i]=x[i]-s*df[i];
	}
	qp_s(xs,Gamma,G_col); // xs=PrX(x-s*df); 
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
	while (fn>f+alpha*(1e-4)*d_aux && k<maxit_r){ // while fn>f+alpha*sig*df'*d && k<=maxit
		printf("    iter_reg=%d fun_n=%f fun=%f\n",k+1,fn,f); //printf ("iter_reg=%d fun_n=%e fun=%e\n",k,fn,f)
		alpha=0.5*alpha; // alpha=ra*alpha;
		for (int i=0;i<G_col;i++){
			xn[i]=x[i]+alpha*d[i]; // xn=x+alpha*d;
		}
		fn=T(xn,df,Delta,c_a,reg_A,G_col,0); // fn=fun(xn);
		k++; // k=k+1;
	} // end
	return k;
}
