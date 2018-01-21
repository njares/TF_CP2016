# Uso:
#
#	genera_D_G(input1,A,p,t_max,r_max)
#
# Donde, "input1" debe ser una direccion de archivo, A debe
# ser una matriz de indices y valores de una matriz de 
# incidencia de un grafo, p y t_max deben ser numeros reales
# y r_max debe ser entero
# input1 = 'g_OD'
# Toma el archivo input1, una matriz de gx2, de pares origen
# destino, y la matriz A de 2ax3, y calcula todas las rutas
# posibles para cada par origen-destino en input1.
# p debe ser el promedio de costo de aristas del grafo dado
# por A, t_max es el tiempo %*máximo*) que %*estará*) buscando nuevas
# rutas %*después*) de haber encontrado la primera, y r_max es la
# cantidad %*máxima*) de rutas que %*encontrará*) por cada par.
# Finalmente, %*generará*) tres archivos:
#   i- 'Delta': un archivo que contiene una matriz de dx2,
#       donde los primero dos elementos de cada fila 
#       describen una %*posición*) en la matriz Delta a donde hay
#       un 1.
#       La matriz asi descrita contiene las rutas posibles
#       entre los pares origen-destino proporcionados por el
#       archivo input1.
#       La matriz no %*contendrá*) ninguna fila para los pares
#       para los que no haya encontrado ninguna ruta posible.
#       El %*tamaño*) de Delta es de tantas filas como rutas haya 
#       encontrado, y tantas columnas como aristas tenga el
#       grafo dado por A.
#       Las rutas %*están*) descriptas con un 1 en cada columna 
#       correspondiente a una arista del grafo que pertenezca
#       a esa ruta.
#   ii- 'Gamma': un archivo que contiene una matriz de gx3,
#       donde los primero dos elementos de cada fila 
#       describen una %*posición*) en la matriz Gamma, y el
#       tercer elemento contiene el valor que %*está*) en esa
#       %*posición*). La matriz asi descrita contiene tantas
#       filas como pares origen-destino a los que se les haya 
#       encontrado al menos una ruta, y tantas columnas como
#       filas tenga Delta. 
#       Cada fila de Gamma contiene un 1 en las columnas que
#       se correspondan con rutas de Delta que sean el par
#       origen-destino correspondiente con la fila de Gamma.
#   iii- 'ignorados': un archivo que contiene una matriz de
#       ix1, de indices de los pares origen-destino a los que
#       no se les haya encontrado ninguna ruta.
#
# NOTA: usa rutas.m
#
function genera_D_G(input1,A,p,t_max,r_max)
	g_OD=load(input1);
	Delta=[];
	Gamma=sparse([],[],[]);
	t=0;
	ignorados=[];
	for i=1:rows(g_OD)
		tic
		printf("Procesando par %d de %d...\n",i,rows(g_OD))
		y=rutas(g_OD(i,1),g_OD(i,2),A,p,t_max,r_max);
		printf("Procesado par %d de %d\n",i,rows(g_OD))
		if (length(y)==0)
			ignorados=[ignorados;i];
			printf("No se encontraron rutas para el par %d, agregado a la lista de pares ignorados\n",i)
			t=t+toc;
			printf("Transcurrieron %d segundos. Tiempo restante estimado: %d segundos\n",t,t*(rows(g_OD)-i)/i)
		else
			i_1=rows(Delta);
			Delta=[Delta;y];
			i_2=rows(Delta);
			t=t+toc;
			printf("Transcurrieron %d segundos. Tiempo restante estimado: %d segundos\n",t,t*(rows(g_OD)-i)/i)
			Gamma(rows(Gamma)+1,rows(Delta))=0;
			Gamma(rows(Gamma),:)=sparse(1,[i_1+1:i_2],1,1,rows(Delta));
		endif
	endfor
	[i j x]=find(Delta);
	dlmwrite("data/output/Delta",[i j x],'	')
	[i j x]=find(Gamma);
	dlmwrite("data/output/Gamma",[i j x],'	')
	dlmwrite("data/output/ignorados",ignorados,'	')
endfunction
