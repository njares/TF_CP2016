#include <stdio.h>

int main(){
	FILE *archivo;
	int X_size=18012,Y_size=3;
	int A[X_size][Y_size];
	char ruta[] = "data/input/matriz_A_2";
	archivo = fopen (ruta, "r");
	
	for (int i=0;i<X_size;i++){
		for (int j=0;j<Y_size;j++){
			fscanf(archivo,"%d",&A[i][j]);
		}
	}
	
	fclose(archivo);
	
	for (int i=0;i<X_size;i++){
		for (int j=0;j<Y_size;j++){
			printf("%d ",A[i][j]);
		}
		printf("\n");
	}
	
	return 0;
}

/*A=load(’data/input/matriz_A_2’);
c_a=load(’data/input/costo_aristas’);
c_a=c_a(find(c_a));
p=mean(c_a);
genera_D_G(’data/input/g_OD’,A,p,3,30);*/
