#include <stdio.h>
#include <stdlib.h>

unsigned int idx(unsigned int i, unsigned int j, unsigned int stride);
void cargar(int * A, char * ruta, int X_size, int Y_size);

int main(){
	unsigned int A_x=18012,A_y=3;
	size_t A_size = A_x * A_y * sizeof(int);
    int * A = malloc(A_size);
	char ruta_A[] = "data/input/matriz_A_2";
	cargar(A,"data/input/matriz_A_2",A_x,A_y); // A=load(’data/input/matriz_A_2’);
	
	for (int i=0;i<A_x;i++){
		for (int j=0;j<A_y;j++){
			printf("%d ",A[idx(i,j,3)]);
		}
		printf("\n");
	}
	return 0;
}

unsigned int idx(unsigned int i, unsigned int j, unsigned int stride) {
    return i * stride + j;
}

void cargar(int * A, char * ruta, int X_size, int Y_size){
	FILE * archivo;
	archivo = fopen (ruta, "r");
	for (int i=0;i<X_size;i++){
		for (int j=0;j<Y_size;j++){
			fscanf(archivo,"%d",&A[idx(i,j,3)]);
		}
	}
	fclose(archivo);
}
	
// c_a=load(’data/input/costo_aristas’);
/*c_a=c_a(find(c_a));
p=mean(c_a);
genera_D_G(’data/input/g_OD’,A,p,3,30);*/
