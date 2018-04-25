#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

unsigned int idx(unsigned int i, unsigned int j, unsigned int stride);
void cargar_int(int * A, char * ruta, unsigned int X_size, unsigned int Y_size);
void cargar_float(float * A, char * ruta, unsigned int X_size, unsigned int Y_size);
void genera_D_G(char * ruta, int * A, float * c_a, float p,float t_max,unsigned int r_max);
unsigned int rutas(int * delta_aux, int i, int j,
	int * A, float * nodos, float p, float t_max, unsigned int r_max, float * c_a);
unsigned int a_star(unsigned int * y_aux, int i, int j, int * A, float * nodos, float * c_a);
unsigned int salientes(int * i_1, int o, int * A);
unsigned int entrantes(int * i_1, int d, int * A);
unsigned int indice_min(float * n_star);
unsigned int arista(unsigned int i_star, unsigned int i_1, int * A);
unsigned int dist_n(int s, int j, float p, float * nodos);

unsigned int A_x=18012,A_y=3,c_a_x=9006,c_a_y=1,n_x=5237,n_y=2;

typedef struct delta {
    int * val;
    unsigned int largo;
    struct delta * next;
} delta_t;


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

void genera_D_G(char * ruta, int * A, float * c_a, float p, float t_max, unsigned int r_max){
	unsigned int g_x=2990,g_y=2,cant_rutas,i_ign=0;
	size_t g_size = g_x * g_y * sizeof(int);
	size_t delta_size = c_a_x * r_max * 2 * sizeof(int);
	size_t nodos_size = n_x * n_y * sizeof(float);
	int * g_OD = malloc(g_size);
	int * delta_aux = malloc(delta_size);
	delta_t * delta_head = NULL; 
	delta_head=malloc(sizeof(delta_t));
	delta_t * current = delta_head; // Delta=[];
	int * ignorados = malloc(g_size); 
	memset(ignorados, -1, g_size); // ignorados=[];
	clock_t end,start;
	float t=0; // t=0;
	cargar_int(g_OD,ruta,g_x,g_y); // g_OD=load(input1);
	float * nodos = malloc(nodos_size);
	cargar_float(nodos,"data/input/nodos_amp",n_x,n_y);	// nodos=load("data/input/nodos_amp");
	for (int i=0;i<g_x;i++){ // for i=1:rows(g_OD)
		start=clock(); // tic
		// printf("Procesando par %d de %d...\n",i,rows(g_OD))
		printf("Procesando par %d de %d...\n",i+1,g_x); 
		// y=rutas(g_OD(i,1),g_OD(i,2),A,p,t_max,r_max);
		memset(delta_aux, 0, delta_size);
		cant_rutas=rutas(delta_aux,g_OD[idx(i,0,g_y)],g_OD[idx(i,1,g_y)],
			A,nodos,p,t_max,r_max,c_a);
		// printf("Procesado par %d de %d\n",i,rows(g_OD))
		printf("Procesado par %d de %d\n",i+1,g_x);
		if (cant_rutas==0){ // if (length(y)==0)
			ignorados[i_ign]=i; // ignorados=[ignorados;i];
			i_ign++;
			// printf("No se encontraron rutas para el par %d, agregado a la lista de pares ignorados\n",i)
			printf("No se encontraron rutas para el par %d, agregado a la lista de pares ignorados\n",i+1);
		} else { // else
			current->val=malloc(delta_size);
			memcpy(current->val,delta_aux,delta_size); // Delta=[Delta;y];
			current->largo=cant_rutas;
			current->val=realloc(current->val,current->largo*sizeof(int));
			current->next=malloc(sizeof(delta_t));
			current->next->next=NULL;
			current=current->next;
		} // endif
		end=clock();
		t+= (float) (end-start)/CLOCKS_PER_SEC; // t=t+toc;
		// printf("Transcurrieron %d segundos. Tiempo restante estimado: %d segundos\n",t,t*(rows(g_OD)-i)/i)
		printf("Transcurrieron %f segundos. Tiempo restante estimado: %f segundos\n",t,t*(g_x-i+1)/(i+1));
	} // endfor
	current=delta_head;
	FILE * delta_file;
	FILE * gamma_file;
	FILE * ignorados_file;
	delta_file=fopen("data/output/Delta","w+");
	gamma_file=fopen("data/output/Gamma","w+");
	ignorados_file=fopen("data/output/ignorados","w+");
/*
	while(current->val != NULL){
		for (int j=0;j<current->largo;j++){
			fprintf(delta_file,"%d ",current->val[j]); // dlmwrite("data/output/Delta",[i j x],'	')
		}
		fprintf(gamma_file,"%d\n",current->largo); // dlmwrite("data/output/Gamma",[i j x],'	')
		fprintf(delta_file,"\n",current->val[j]);
		current=current->next;
	}
*/
	for (int i=0;i<i_ign;i++){
		fprintf(ignorados_file,"%d\n",ignorados[i]); // dlmwrite("data/output/ignorados",ignorados,'	')
	}
	while(delta_head->next != NULL){
		current=delta_head;
		delta_head=delta_head->next;
		free(current->val);
		free(current);
	}
	free(g_OD);
	free(ignorados);
	free(delta_aux);
} // endfunction

unsigned int rutas(int * delta_aux, int i,int j,
	int * A, float * nodos, float p, float t_max,unsigned int r_max, float * c_a)
{
	unsigned int largo_y_aux,n_max,y_aux[100],i_n=0,i_1[30],i_1_size,i_s=0,current,i_aux;
	float t;
	largo_y_aux=a_star(y_aux,i,j,A,nodos,c_a); // y_aux=a_star(i,j,A); % busca el camino mas corto con A*
	if (largo_y_aux==0){ // if (length(y_aux)==0)
		return 0; // y=[];
	} else { // else
		// % calcula la distancia maxima de camino
		n_max=(int)(((float)largo_y_aux)*1.25)+1; // n_max=floor(length(y_aux)*1.25)+1; 
		// y_aux(n_max)=0;
		size_t size_n=n_max*r_max*sizeof(unsigned int);
		size_t size_s=n_max*n_max*30*sizeof(unsigned int);
		int * n=malloc(size_n); // n=[];
		int * s=malloc(size_s); // s=[];
		memset(n,-1,size_n);
		memset(s,-1,size_s);
		memcpy(n,y_aux,largo_y_aux*sizeof(unsigned int)); // n=[n;y_aux];
		i_n++;
		if (i>0){ // if(i>0)
			// % defino la pila para dfs como si el primer camino encontrado
			// % hubiera sido el encontrado con a_star
			// k_y=length(find(y_aux));
			for (int k=0;k<largo_y_aux-1;k++) { // for k=1:k_y-1
				// % encuentro salientes para nodos salvo el %*último*)
				i_1_size=salientes(i_1,y_aux[k],A); // i_1=salientes(y_aux(k_y-k),A);
				for (int l=0;l<i_1_size;l++) {
					// [c i1 i2]=intersect(i_1,y_aux(k_y-k+1));
					// % saco el que ya %*había*) sido 'encontrado'
					if (i_1[l]!=y_aux[k+1]) { // i_1=i_1(setdiff(1:rows(i_1),i2)); 
						// % calculo los caminos para agregar a la pila
						for (int m=0;m<k+1;m++) { 
							s[idx(i_s,m,n_max)]=y_aux[m]; // v=y_aux(1:(k_y-k));
							// v(rows(v),length(y_aux))=0;
							// v=repmat(v,length(i_1),1);
						}
						s[idx(i_s,k+1,n_max)]=i_1[l]; // // v(:,k_y-k+1)=i_1; 
						// % agrego los nuevos caminos a la pila
						i_s++; // s=[v;s]; 
					}
				}
			} // endfor
/*
			k=0;
			l=0;
			i_max=0;
			i_min=n_max+1;
			tic
*/
			while (i_s!=0 && t<t_max && i_n < r_max) { // while (rows(s)!=0 && toc<t_max && rows(n)<r_max)
				// % pongo el ultimo elemento de la pila en v
				current=i_s; // v=s(rows(s),:);	
				// % saco el ultimo elemento de la pila
				i_s--; // s=s([1:rows(s)-1],:);
				// % guardo la posicion del %*último*) elemento de v
				i_aux=-1;
				for (k=0;k<n_max;k++) {
					if (s[idx(current,k,n_max)]==-1) {
						i_aux=k-1; // i_aux=max(find(v));
						k=n_max;
					}
				}
				if (i_aux==-1) {
					i_aux=n_max-1;
				}
				if (fabs(nodos[idx(s[idx(current,i_aux,n_max)],0,n_y)]-nodos[idx(j-1,0,n_y)])*94870 +
					fabs(nodos[idx(s[idx(current,i_aux,n_max)],1,n_y)]-nodos[idx(j-1,1,n_y)])*111180 < 300) { // if (sum(abs(nodos(v(i_aux),:)-nodos(j,:)).*[94870 111180])<300) 
					// % si estoy a menos de 300 m del destino, agrego el camino a n
					memcpy(&n[idx(i_n,0,n_max)],&s[idx(current,0,n_max)],n_max*sizeof(int)); // n=[n;v];
/*
					k=k+1;
					i_max=max(i_max,i_aux);
					i_min=min(i_min,i_aux);
*/
				} else if (i_aux==n_max-1) { // elseif (i_aux==n_max)
					// % si no llegue, y alcance la cantidad %*máxima*) de nodos, sigo 
					// % con el proximo
/*
					k=k+1;
*/					
					continue; // continue;
				} else if (n_max-i_aux<dist_n(s[idx(current,i_aux,n_max)],j,p,nodos)) { // elseif (n_max-i_aux<dist_n(v(i_aux),j,p,nodos))
					// % si no llegue, y me %*alejé*) demasiado, descarto este camino
/*
					l=l+2**(n_max-i_aux);
*/
					continue; // continue;
				} else { // else
					// % si no llegue, y no alcance la cantidad %*máxima*) de nodos, 
					// % agrego nuevos nodos
					// % averiguo los adyacentes al %*último*) elemento de v
					i_1_size=salientes(i_1,s[idx(current,i_aux,n_max)],A); // i_1=salientes(v(i_aux),A);
					// % chequeo si se repite con alguno de los ya visitados
					int flag;
					for (int k=0;k<i_1_size;k++) {
						// [c i1 i2]=intersect(v,i_1);
						// % saco los repetidos
						flag=0;
						for (int l=0;l<i_aux && !flag;l++) {
							if (i_1[k]==s[idx(current,l,n_max)]) { // i_1=i_1(setdiff(1:rows(i_1),i2)); 
								flag=1;
							}
						}
						if (flag) {
							continue;
						} else {
							// % calculo los nuevos caminos para agregar a la pila
							for (int l=0;l<i_aux+1;l++) {
								s[idx(i_s,l,n_max)]=s[idx(current,l,n_max)];
								// v=repmat(v,length(i_1),1);
							}
							s[idx(i_s,i_aux+1,n_max)]=i_1[k]; // // v(:,i_aux+1)=i_1;
							// % agrego los nuevos caminos a la pila
							i_s++; // s=[s;v];
						}
					}
				} // endif
			} // endwhile
		}else{ // else % Lo mismo, pero 'contramano'
			i=-i; // i=-i;
			j=-j; // j=-j;
			// k_y=length(find(y_aux));
			for (int k=largo_y_aux-1;k>0;k--) { // for k=2:k_y
				i_1_size=entrantes(i_1,y_aux[k],A); // i_1=entrantes(y_aux(k),A);
				for (int l=0;l<i_1_size;l++) {
					// [c i1 i2]=intersect(i_1,y_aux(k-1));
					if (i_1[l]!=y_aux[k-1]) { // i_1=i_1(setdiff(1:rows(i_1),i2));
						for (int m=1;m<largo_y_aux-k+1;m++) {
							s[idx(i_s,m,n_max)]=y_aux[k+m-1]; // v=y_aux(k:k_y);
							// v(rows(v),length(y_aux)-1)=0;
							// v=[0,v];
							// v=repmat(v,length(i_1),1);
						}
						s[idx(i_s,0,n_max)]=i_1[l]; // v(:,1)=i_1;
						i_s++; // s=[v;s];
					}
				}
			} // endfor 
/*
			k=0;
			l=0;
			i_max=0;
			i_min=n_max+1;
			tic
*/
			while (i_s!=0 && t<t_max && i_n < r_max) { // while (rows(s)!=0 && toc<t_max && rows(n)<r_max)
				current=i_s; // v=s(rows(s),:);
				i_s--; // s=s([1:rows(s)-1],:);
				i_aux=-1;
				for (k=0;k<n_max;k++) {
					if (s[idx(current,k,n_max)]==-1) {
						i_aux=k-1; // i_aux=max(find(v));
						k=n_max;
					}
				}
				if (i_aux==-1) {
					i_aux=n_max-1;
				}
				if (fabs(nodos[idx(s[idx(current,0,n_max)],0,n_y)]-nodos[idx(i-1,0,n_y)])*94870 +
					fabs(nodos[idx(s[idx(current,0,n_max)],1,n_y)]-nodos[idx(i-1,1,n_y)])*111180 < 300) { // if (sum(abs(nodos(v(1),:)-nodos(i,:)).*[94870 111180])<300)
					memcpy(&n[idx(i_n,0,n_max)],&s[idx(current,0,n_max)],n_max*sizeof(int)); // n=[n;v];
/*
					k=k+1;
					i_max=max(i_max,i_aux);
					i_min=min(i_min,i_aux);
*/
				} else if (i_aux==n_max-1) { // elseif (i_aux==n_max)
/*
					k=k+1;
*/					
					continue; // continue;
				} else if (n_max-i_aux<dist_n(s[idx(current,0,n_max)],i,p,nodos)) { // elseif (n_max-i_aux<dist_n(v(1),i,p,nodos))
/*
					l=l+2**(n_max-i_aux);
*/
					continue; // continue;
				} else { // else
					i_1_size=entrantes(i_1,s[idx(current,0,n_max)],A); // i_1=entrantes(v(1),A); 
					int flag;
					for (int k=0;k<i_1_size;k++) {
						// [c i1 i2]=intersect(v,i_1);
						flag=0;
						for (int l=0;l<i_aux && !flag;l++) {
							if (i_1[k]==s[idx(current,l,n_max)]) { // i_1=i_1(setdiff(1:rows(i_1),i2)); 
								flag=1;
							}
						}
						if (flag) {
							continue;
						} else {
							for (int l=1;l<i_aux+2;l++){
								s[idx(i_s,l,n_max)]=s[idx(current,l-1,n_max)];
								// v=[0 v(1:(length(v)-1))];
								// v=repmat(v,length(i_1),1);
							}
							s[i_s,0;n_max]=i_1[k]; // v(:,1)=i_1;
							i_s++; // s=[s;v];
						}
					}
				} // endif
			} // endwhile			
		} // endif

/*
		y=[];
% reconstruyo los caminos para devolverlos
		for k=1:rows(n)
			y_aux=sparse([],[],[],1,max(A(:,1))); 
			for l=1:max(find(n(k,:)))-1
				y_aux=y_aux+arista(n(k,l),n(k,l+1),A);
			endfor
			y=[y;y_aux];
		endfor
*/

		return 1;
	} // endif
	free(n);
	free(s);
} // endfunction

unsigned int a_star(unsigned int * y_aux, int o, int d, int * A, float * nodos, float * c_a){
	size_t n_size = n_x * 4 * sizeof(float);
	unsigned int i_1[30],i_1_size,open_size=0,i_star,flag,ruta_reves[100],i_r=0;
	float * n_star = malloc(n_size); // [conjunto costo_acum costo_falt came_from]
	float h_aux,current[4],g_aux;
	for (int i=0;i<n_x;i++) {
		n_star[idx(i,0,4)]=-1;
		n_star[idx(i,1,4)]=-1;
		n_star[idx(i,2,4)]=-1;
		n_star[idx(i,3,4)]=-1;
	} // closed=[]; // came_from=zeros(rows(nodos),1);
	if (o>0) { // if(o>0) % si el grafo debe recorrerse 'mano'
		// % si origen y destino %*están*) a menos de 300 m
		if ( fabs(nodos[idx(o-1,0,n_y)]-nodos[idx(d-1,0,n_y)])*94870 +
			 fabs(nodos[idx(o-1,1,n_y)]-nodos[idx(d-1,1,n_y)])*111180 < 300) { // if (sum(abs(nodos(o,:)-nodos(d,:)).*[94870 111180])<300) 
			i_1_size=salientes(i_1,o,A); // i_1=(salientes(o,A))(1);
			y_aux[0]=o;
			y_aux[1]=i_1[0]; // y=[o i_1];
			return 2;
		} else { // else
			h_aux=fabs(nodos[idx(o-1,0,n_y)]-nodos[idx(d-1,0,n_y)]) +
					fabs(nodos[idx(o-1,1,n_y)]-nodos[idx(d-1,1,n_y)]); // h_aux=sum(abs(nodos(o,:)-nodos(d,:)));
			n_star[idx(o-1,0,4)]=1;
			n_star[idx(o-1,1,4)]=0;
			n_star[idx(o-1,2,4)]=h_aux;
			open_size++; // open=[o 0 h_aux];
			while (open_size!=0) { // while (rows(open)!=0)
				i_star=indice_min(n_star); // [x i]=min(open(:,2)+open(:,3));
				// % guardo el nodo mas barato en current
				memcpy(current,&n_star[idx(i_star,0,4)],4*sizeof(float)); // current=open(i,:);
				if ( fabs(nodos[idx(i_star,0,n_y)]-nodos[idx(d-1,0,n_y)])*94870 +
					 fabs(nodos[idx(i_star,1,n_y)]-nodos[idx(d-1,1,n_y)])*111180 < 300) { // if (sum(abs(nodos(current(1),:)-nodos(d,:)).*[94870 111180])<300) 
					// % si estoy a menos de 300 m del destino, termino
					break; // break;
				} // endif
				// % agrego el nodo actual a closed
				// closed=[closed;open(i,:)];
				// % saco el nodo actual de open
				// i_aux=setdiff([1:rows(open)],i);
				n_star[idx(i_star,0,4)]=2; 
				open_size--; // open=open(i_aux,:);
				// % busco los vecinos de current
				i_1_size=salientes(i_1,i_star+1,A); // i_1=salientes(current(1),A);
				for (int i=0;i<i_1_size;i++) { // for i=1:rows(i_1)
					// % si ya esta en closed, sigo
					if (n_star[idx(i_1[i]-1,0,4)] == 2) { // if (length(find(closed(:,1)==i_1(i))))
						continue; // continue;
					} // endif
					g_aux=current[1]+c_a[arista(i_star,i_1[i],A)-1]; // g_aux=current(2)+c_a(find(arista(current(1),i_1(i),A)));
					h_aux=fabs(nodos[idx(i_1[i]-1,0,n_y)]-nodos[idx(d-1,0,n_y)]) +
						  fabs(nodos[idx(i_1[i]-1,1,n_y)]-nodos[idx(d-1,1,n_y)]); // h_aux=sum(abs(nodos(i_1(i),:)-nodos(d,:)));
					if (n_star[idx(i_1[i]-1,0,4)] == 1) { // if (length(find(open(:,1)==i_1(i)))) 
						// % si ya %*está*) en open, lo dejo con el menor costo
						// i_aux=find(open(:,1)==i_1(i));
						if (n_star[idx(i_1[i]-1,1,4)]>g_aux) { // if(open(i_aux,2)>g_aux)
							n_star[idx(i_1[i]-1,1,4)]=g_aux; // open(i_aux,2)=g_aux;
							n_star[idx(i_1[i]-1,3,4)]=i_star; // came_from(open(i_aux,1))=current(1);
						} // endif
					} else { // else
						// % si no %*está*) en open, lo agrego
						open_size++;
						n_star[idx(i_1[i]-1,0,4)]=1;
						n_star[idx(i_1[i]-1,1,4)]=g_aux;
						n_star[idx(i_1[i]-1,2,4)]=h_aux; // open=[open;[i_1(i) g_aux h_aux]];
						n_star[idx(i_1[i]-1,3,4)]=i_star; // came_from(i_1(i))=current(1);
					} // endif
				} // endfor
			} // endwhile
			if (open_size==0) { // if (rows(open)==0)
				return 0; // y=[];
			} else { // else
				// % reconstruyo el camino
				flag=i_star; // flag=current(1);
				ruta_reves[i_r]=i_star; // y=current(1);
				i_r++;
				while (flag!=o-1) { // while (flag!=o)
					ruta_reves[i_r]=n_star[idx(flag,3,4)];
					i_r++; // y=[came_from(flag);y];
					flag=n_star[idx(flag,3,4)]; // flag=came_from(flag);
				} // endwhile
				// y=y';
				for (int i=0;i<i_r;i++) {
					y_aux[i]=ruta_reves[i_r-i-1];
				}
				return i_r;
			} // endif
		} // endif
	} else { // else % si el grafo debe recorrerse 'contramano'
		o=-o; // o=-o;
		d=-d; // d=-d;
		if ( fabs(nodos[idx(o-1,0,2)]-nodos[idx(d-1,0,2)])*94870 +
			 fabs(nodos[idx(o-1,1,2)]-nodos[idx(d-1,1,2)])*111180 < 300) { // if (sum(abs(nodos(o,:)-nodos(d,:)).*[94870 111180])<300) 
			i_1_size=entrantes(i_1,d,A); // i_1=(entrantes(d,A))(1);
			y_aux[0]=i_1[0];
			y_aux[1]=d; // y=[i_1 d];
			return 2;
		} else { // else
			h_aux=fabs(nodos[idx(o-1,0,n_y)]-nodos[idx(d-1,0,n_y)]) +
					fabs(nodos[idx(o-1,1,n_y)]-nodos[idx(d-1,1,n_y)]); // h_aux=sum(abs(nodos(o,:)-nodos(d,:)));
			n_star[idx(d-1,0,4)]=1;
			n_star[idx(d-1,1,4)]=0;
			n_star[idx(d-1,2,4)]=h_aux;
			open_size++; // open=[d 0 h_aux];
			while (open_size!=0) { // while (rows(open)!=0)
				i_star=indice_min(n_star); // [x i]=min(open(:,2)+open(:,3));
				memcpy(current,&n_star[idx(i_star,0,4)],4*sizeof(float)); // current=open(i,:);
				if ( fabs(nodos[idx(i_star,0,n_y)]-nodos[idx(o-1,0,n_y)])*94870 +
					 fabs(nodos[idx(i_star,1,n_y)]-nodos[idx(o-1,1,n_y)])*111180 < 300) { // if (sum(abs(nodos(current(1),:)-nodos(o,:)).*[94870 111180])<300) 
					break; // break;
				} // endif
				// closed=[closed;open(i,:)];
				// i_aux=setdiff([1:rows(open)],i);
				n_star[idx(i_star,0,4)]=2; 
				open_size--; // open=open(i_aux,:);
				i_1_size=entrantes(i_1,i_star+1,A); // i_1=entrantes(current(1),A);
				for (int i=0;i<i_1_size;i++) { // for i=1:rows(i_1)
					if (n_star[idx(i_1[i]-1,0,4)] == 2) { // if (length(find(closed(:,1)==i_1(i))))
						continue; // continue;
					} // endif
					g_aux=current[1]+c_a[arista(i_star,i_1[i],A)-1]; // g_aux=current(2)+c_a(find(arista(i_1(i),current(1),A)));
					h_aux=fabs(nodos[idx(i_1[i]-1,0,n_y)]-nodos[idx(o-1,0,n_y)]) +
						  fabs(nodos[idx(i_1[i]-1,1,n_y)]-nodos[idx(o-1,1,n_y)]); // h_aux=sum(abs(nodos(i_1(i),:)-nodos(o,:)));
					if (n_star[idx(i_1[i]-1,0,4)] == 1) { // if (length(find(open(:,1)==i_1(i)))) 
						if (n_star[idx(i_1[i]-1,1,4)]>g_aux) { // if(open(i_aux,2)>g_aux)
							n_star[idx(i_1[i]-1,1,4)]=g_aux; // open(i_aux,2)=g_aux;
							n_star[idx(i_1[i]-1,3,4)]=i_star; // came_from(open(i_aux,1))=current(1);
						} // endif
					} else { // else
						open_size++;
						n_star[idx(i_1[i]-1,0,4)]=1;
						n_star[idx(i_1[i]-1,1,4)]=g_aux;
						n_star[idx(i_1[i]-1,2,4)]=h_aux; // open=[open;[i_1(i) g_aux h_aux]];
						n_star[idx(i_1[i]-1,3,4)]=i_star; // came_from(i_1(i))=current(1);
					} // endif
				} // endfor
			} // endwhile
			if (open_size==0) { // if (rows(open)==0)
				return 0; // y=[];
			} else { // else
				flag=i_star; // flag=current(1);
				y_aux[i_r]=i_star; // y=current(1);
				i_r++;
				while (flag!=d-1) { // while (flag!=d)
					y_aux[i_r]=n_star[idx(flag,3,4)];
					i_r++; // y=[y;came_from(flag)];
					flag=n_star[idx(flag,3,4)]; // flag=came_from(flag);
				} // endwhile
				// y=y';
				return i_r;
			} // endif
		} // endif
	} // endif
	free(n_star);
} // endfunction

unsigned int salientes(int * i_1, int o, int * A){
	int i_aux[30],size=0,k=0;
	for (int i=0;i<A_x;i++){
		if (A[idx(i,1,A_y)]==o &&
			A[idx(i,2,A_y)]==1 ){
			i_aux[size]=A[idx(i,0,A_y)];
			size++;
		}
	}
	for (int i=0;i<A_x;i++){
		for (int j=0;j<size;j++){
			if ( A[idx(i,0,A_y)]==i_aux[j] && 
				 A[idx(i,1,A_y)]!=o) {
				i_1[k]=A[idx(i,1,A_y)];
				k++;
			}
		}
	}
	return size;
}

unsigned int entrantes(int * i_1, int d, int * A){
	int i_aux[30],size=0,k=0;
	for (int i=0;i<A_x;i++){
		if (A[idx(i,1,A_y)]==d &&
			A[idx(i,2,A_y)]==-1 ){
			i_aux[size]=A[idx(i,0,A_y)];
			size++;
		}
	}
	for (int i=0;i<A_x;i++){
		for (int j=0;j<size;j++){
			if ( A[idx(i,0,A_y)]==i_aux[j] && 
				 A[idx(i,1,A_y)]!=d) {
				i_1[k]=A[idx(i,1,A_y)];
				k++;
			}
		}
	}
	return size;
}

unsigned int indice_min(float * n_star){
	unsigned int i_min;
	float min=100000000;
	for (int i=0;i<n_x;i++){
		if (n_star[idx(i,0,4)]==1) {
			if (min>n_star[idx(i,1,4)]+n_star[idx(i,2,4)]) {
				min=n_star[idx(i,1,4)]+n_star[idx(i,2,4)];
				i_min=i;
			}
		}
	}
	return i_min;
}

unsigned int arista(unsigned int i_star, unsigned int i_1, int * A){
	unsigned int a_1[100],a_2[100],i_a_1=0,i_a_2=0;
	for (int i=0;i<A_x;i++){
		if (A[idx(i,1,A_y)]==i_star+1) {
			a_1[i_a_1]=A[idx(i,0,A_y)];
			i_a_1++;
		}
		if (A[idx(i,1,A_y)]==i_1) {
			a_2[i_a_2]=A[idx(i,0,A_y)];
			i_a_2++;
		}
	}
	for (int i=0;i<i_a_1;i++){
		for (int j=0;j<i_a_2;j++) {
			if (a_1[i]==a_2[j]) {
				return a_1[i];
			}
		}
	}
	return 0;
}

unsigned int dist_n(int s, int j, float p, float * nodos){
	float dist;
	dist=fabs(nodos[idx(s-1,0,n_y)]-nodos[idx(j-1,0,n_y)]) +
		 fabs(nodos[idx(s-1,1,n_y)]-nodos[idx(j-1,1,n_y)]); // dist=sum(abs(nodos(i,:)-nodos(j,:)));
	return (int)(dist/p) + 4; // y=floor(dist/p)+4;
}
