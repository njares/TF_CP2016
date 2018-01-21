#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>

unsigned int idx(unsigned int i, unsigned int j, unsigned int stride);
void cargar_int(int * A, char * ruta, unsigned int X_size, unsigned int Y_size);
void cargar_float(float * A, char * ruta, unsigned int X_size, unsigned int Y_size);
void genera_D_G(char * ruta, int * A, float * c_a, float p,float t_max,unsigned int r_max);
unsigned int rutas(int * Delta_aux, int i,int j,int * A, float p,float t_max,unsigned int r_max);

unsigned int A_x=18012,A_y=3,c_a_x=9006,c_a_y=1;

int main(){
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
	unsigned int r_max=30;
	float t_max=3;
	//genera_D_G(’data/input/g_OD’,A,p,3,30);
	genera_D_G("data/input/g_OD",A,c_a,p,t_max,r_max);

	free(A);
	free(c_a);
	return 0;
}

unsigned int idx(unsigned int x, unsigned int y, unsigned int stride) {
    return x * stride + y;
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

void genera_D_G(char * ruta, int * A,float * c_a,float p,float t_max,unsigned int r_max){
	unsigned int g_x=2990,g_y=2,cant_rutas,i_ign=0;
	size_t g_size = g_x * g_y * sizeof(int);
	size_t delta_size = c_a_x * r_max * 2 * sizeof(int);
	int * g_OD = malloc(g_size);
	int * delta_aux = malloc(delta_size); // Delta=[];
// usar linked lists para los arreglos incrementales, total no tengo que 
// accederlos nunca salvo para escribirlos al final	
	int * ignorados = malloc(g_size); // ignorados=[];
	clock_t end,start;
	float t=0; // t=0;
	memset(ignorados, -1, g_size);
	memset(delta_aux, 0, delta_size);
	// g_OD=load(input1);
	cargar_int(g_OD,ruta,g_x,g_y);
	for (int i=0;i<g_x;i++){ // for i=1:rows(g_OD)
		start=clock(); // tic
		// printf("Procesando par %d de %d...\n",i,rows(g_OD))
		printf("Procesando par %d de %d...\n",i+1,g_x); 
		// y=rutas(g_OD(i,1),g_OD(i,2),A,p,t_max,r_max);
		cant_rutas=rutas(delta_aux,g_OD[idx(i,0,g_y)],g_OD[idx(i,1,g_y)],A,p,t_max,r_max);
		// printf("Procesado par %d de %d\n",i,rows(g_OD))
		printf("Procesado par %d de %d\n",i+1,g_x);
		if (cant_rutas==0){ // if (length(y)==0)
			ignorados[i_ign]=i; // ignorados=[ignorados;i];
			i_ign++;
			// printf("No se encontraron rutas para el par %d, agregado a la lista de pares ignorados\n",i)
			printf("No se encontraron rutas para el par %d, agregado a la lista de pares ignorados\n",i+1);
		} else { // else
/*
			Delta=[Delta;y];
*/
		} // endif
		end=clock();
		t+= (float) (end-start)/CLOCKS_PER_SEC; // t=t+toc;
		// printf("Transcurrieron %d segundos. Tiempo restante estimado: %d segundos\n",t,t*(rows(g_OD)-i)/i)
		printf("Transcurrieron %f segundos. Tiempo restante estimado: %f segundos\n",t,t*(g_x-i+1)/(i+1));
	} // endfor
/*
			i_1=rows(Delta);
			i_2=rows(Delta);
			Gamma(rows(Gamma)+1,rows(Delta))=0;
			Gamma(rows(Gamma),:)=sparse(1,[i_1+1:i_2],1,1,rows(Delta));	
	Gamma=sparse([],[],[]);
	[i j x]=find(Delta);
	dlmwrite("data/output/Delta",[i j x],'	')
	[i j x]=find(Gamma);
	dlmwrite("data/output/Gamma",[i j x],'	')
	dlmwrite("data/output/ignorados",ignorados,'	')
*/	
	free(g_OD);
	free(ignorados);
	free(delta_aux);
} // endfunction

unsigned int rutas(int * Delta_aux, int i,int j,int * A, float p,float t_max,unsigned int r_max){
	sleep(1);
	return 0;
}

/*
# Uso:
#
#	y=rutas(i,j,A,p,t_max,r_max)
#
# Donde i, j y r_max deben ser enteros, A debe ser una matriz
# de indices y valores de una matriz de incidencia de un 
# grafo, y p y t_max deben ser reales.
# Dados i y j indices de nodos de un grafo dado por la matriz
# A de nxm, devuelve el arreglo y, de rxn, donde r es la 
# cantidad de rutas encontradas entre i y j, y n es la 
# cantidad de aristas del grafo. 'p' es el promedio de costo
# de aristas del grafo, la %*ejecución*) de la segunda parte del 
# algoritmo no excede t_max, y no busca mas de r_max rutas.
# Primero realiza una %*búsqueda*) con un algoritmo 'A estrella',
# y luego lleva adelante un algoritmo 'Busqueda en 
# profundidad' para encontrar %*más*) rutas.
# La cantidad %*máxima*) de nodos que puede tener un camino esta
# en %*función*) de la longitud del camino %*óptimo*) encontrado con
# 'A estrella' (se permiten aproximadamente un 25 % de nodos
# %*más*)).
# El promedio de costo de aristas de A es usado al llamar 
# a dist_n para descartar caminos que se hayan alejado mucho
# del destino.
# Las rutas (filas) de y, son un arreglo que tiene un 1 en la
# columna que se corresponde con cada arista que haya sido 
# usada en esa ruta.
# 'y' %*está*) en formato sparse
# Si no encuentra %*ningún*) camino, devuelve un vector %*vacío*).
#
# NOTA: usa a_star.m, salientes.m, entrantes.m, dist_n.m y 
# arista.m
#
function y=rutas(i,j,A,p,t_max,r_max)
	nodos=load("data/input/nodos_amp");
	y_aux=a_star(i,j,A); % busca el camino mas corto con A*
	if (length(y_aux)==0)
		y=[];
	else
		n=[];
% calcula la distancia maxima de camino
		n_max=floor(length(y_aux)*1.25)+1; 
		y_aux(n_max)=0;
		n=[n;y_aux];
		s=[];
		if(i>0)
% defino la pila para dfs como si el primer camino encontrado
% hubiera sido el encontrado con a_star
			k_y=length(find(y_aux));
			for k=1:k_y-1
% encuentro salientes para nodos salvo el %*último*)
				i_1=salientes(y_aux(k_y-k),A); 
				[c i1 i2]=intersect(i_1,y_aux(k_y-k+1));
% saco el que ya %*había*) sido 'encontrado'
				i_1=i_1(setdiff(1:rows(i_1),i2)); 
				v=y_aux(1:(k_y-k));
				v(rows(v),length(y_aux))=0;
				v=repmat(v,length(i_1),1);
% calculo los caminos para agregar a la pila
				v(:,k_y-k+1)=i_1; 
% agrego los nuevos caminos a la pila
				s=[v;s]; 
			endfor
			k=0;
			l=0;
			i_max=0;
			i_min=n_max+1;
			tic
			while (rows(s)!=0 && toc<t_max && rows(n)<r_max)
% pongo el ultimo elemento de la pila en v
				v=s(rows(s),:);	
% saco el ultimo elemento de la pila
				s=s([1:rows(s)-1],:);	
% guardo la posicion del %*último*) elemento de v
				i_aux=max(find(v));		
				if (sum(abs(nodos(v(i_aux),:)-nodos(j,:)).*[94870 111180])<300) 
% si estoy a menos de 300 m del destino, agrego el camino a n
					n=[n;v];
					k=k+1;
					i_max=max(i_max,i_aux);
					i_min=min(i_min,i_aux);
				elseif (i_aux==n_max)
% si no llegue, y alcance la cantidad %*máxima*) de nodos, sigo 
% con el proximo
					k=k+1;
					continue;
				elseif (n_max-i_aux<dist_n(v(i_aux),j,p,nodos))
% si no llegue, y me %*alejé*) demasiado, descarto este camino
					l=l+2**(n_max-i_aux);
					continue;
				else
% si no llegue, y no alcance la cantidad %*máxima*) de nodos, 
% agrego nuevos nodos
% averiguo los adyacentes al %*último*) elemento de v
					i_1=salientes(v(i_aux),A);	
% chequeo si se repite con alguno de los ya visitados
					[c i1 i2]=intersect(v,i_1);	
% saco los repetidos
					i_1=i_1(setdiff(1:rows(i_1),i2)); 
					v=repmat(v,length(i_1),1);
% calculo los nuevos caminos para agregar a la pila
					v(:,i_aux+1)=i_1;	
% agrego los nuevos caminos a la pila			
					s=[s;v];			
				endif
			endwhile
		else % Lo mismo, pero 'contramano'
			i=-i;
			j=-j;
			k_y=length(find(y_aux));
			for k=2:k_y
				i_1=entrantes(y_aux(k),A);
				[c i1 i2]=intersect(i_1,y_aux(k-1));
				i_1=i_1(setdiff(1:rows(i_1),i2)); 
				v=y_aux(k:k_y);
				v(rows(v),length(y_aux)-1)=0;
				v=[0,v];
				v=repmat(v,length(i_1),1);
				v(:,1)=i_1;
				s=[v;s];
			endfor
			k=0;
			l=0;
			i_max=0;
			i_min=n_max+1;
			tic
			while (rows(s)!=0 && toc<t_max && rows(n)<r_max)
				v=s(rows(s),:);	
				s=s([1:rows(s)-1],:);
				i_aux=max(find(v));
				if (sum(abs(nodos(v(1),:)-nodos(i,:)).*[94870 111180])<300) 
					n=[n;v];
					k=k+1;
					i_max=max(i_max,i_aux);
					i_min=min(i_min,i_aux);
				elseif (i_aux==n_max)
					k=k+1;
					continue;
				elseif (n_max-i_aux<dist_n(v(1),i,p,nodos))
					l=l+2**(n_max-i_aux);
					continue;
				else
					i_1=entrantes(v(1),A); 
					[c i1 i2]=intersect(v,i_1);
					i_1=i_1(setdiff(1:rows(i_1),i2));
					v=[0 v(1:(length(v)-1))];
					v=repmat(v,length(i_1),1);
					v(:,1)=i_1;
					s=[s;v];
				endif
			endwhile
		endif
		y=[];
% reconstruyo los caminos para devolverlos
		for k=1:rows(n)
			y_aux=sparse([],[],[],1,max(A(:,1))); 
			for l=1:max(find(n(k,:)))-1
				y_aux=y_aux+arista(n(k,l),n(k,l+1),A);
			endfor
			y=[y;y_aux];
		endfor
	endif
endfunction
*/
