#include <stdio.h>
#include <stdlib.h>

unsigned int idx(unsigned int i, unsigned int j, unsigned int stride);
void cargar_int(int * A, char * ruta, int X_size, int Y_size);
void cargar_float(float * A, char * ruta, int X_size, int Y_size);

int main(){
	unsigned int A_x=18012,A_y=3,c_a_x=9006,c_a_y=1;
	size_t A_size = A_x * A_y * sizeof(int);
	size_t c_a_size = c_a_x * c_a_y * sizeof(float);
    int * A = malloc(A_size);
    float * c_a = malloc(c_a_size);
    // A=load(’data/input/matriz_A_2’);
	cargar_int(A,"data/input/matriz_A_2",A_x,A_y);
	// c_a=load(’data/input/costo_aristas’);	
	cargar_float(c_a,"data/input/costo_aristas",c_a_x,c_a_y);
	// c_a=c_a(find(c_a));
	// p=mean(c_a);
	float p=0;
	int n=0;
	for (int i=0;i<c_a_x;i++){
		for (int j=0;j<c_a_y;j++){
			if (c_a[idx(i,j,c_a_y)]!=0.0){
				p+=c_a[idx(i,j,c_a_y)];
				n++;
			}
		}
	}
	p=p/n;
	printf ("%f\n",p);
	
/*	for (int i=0;i<A_x;i++){
		for (int j=0;j<A_y;j++){
			printf("%d ",A[idx(i,j,3)]);
		}
		printf("\n");
	}
	getchar();
	for (int i=0;i<c_a_x;i++){
		for (int j=0;j<c_a_y;j++){
			printf("%d ",c_a[idx(i,j,3)]);
		}
		printf("\n");
	}*/
	return 0;
}

unsigned int idx(unsigned int x, unsigned int y, unsigned int stride) {
    return x * stride + y;
}

void cargar_int(int * A, char * ruta, int X_size, int Y_size){
	FILE * archivo;
	archivo = fopen (ruta, "r");
	for (int i=0;i<X_size;i++){
		for (int j=0;j<Y_size;j++){
			fscanf(archivo,"%d",&A[idx(i,j,Y_size)]);
		}
	}
	fclose(archivo);
}

void cargar_float(float * A, char * ruta, int X_size, int Y_size){
	FILE * archivo;
	archivo = fopen (ruta, "r");
	for (int i=0;i<X_size;i++){
		for (int j=0;j<Y_size;j++){
			fscanf(archivo,"%f",&A[idx(i,j,Y_size)]);
		}
	}
	fclose(archivo);
}

//genera_D_G(’data/input/g_OD’,A,p,3,30);
