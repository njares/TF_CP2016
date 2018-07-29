#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

void cargar_int(int * A, char * ruta, unsigned int X_size, unsigned int Y_size);
void cargar_float(float * A, char * ruta, unsigned int X_size, unsigned int Y_size);
unsigned int idx(unsigned int i, unsigned int j, unsigned int stride);

unsigned int c_a_x=9006,c_a_y=1,D_x=355746,D_y=2,G_x=2990,G_y=1;

int main(){
	size_t c_a_size = c_a_x * c_a_y * sizeof(float);
	size_t D_size = D_x * D_y * sizeof(int);
	size_t G_size = G_x * G_y * sizeof(int);
	int reg_A[4];
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
