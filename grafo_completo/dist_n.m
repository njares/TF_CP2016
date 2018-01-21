# Uso:
#
#	y=dist_n(i,j,p,nodos)
#
# Donde i y j deben ser enteros, p debe ser real, y nodos
# debe ser una matriz de nx2.
# Dados i y j indices de nodos de un grafo, p el promedio de
# costo de aristas, y nodos la matriz con todas las 
# coordenadas de los nodos de un grafo, devuelve una cota 
# %*máxima*) de cantidad de aristas para ir de i a j. 
#
function y=dist_n(i,j,p,nodos)
	dist=sum(abs(nodos(i,:)-nodos(j,:)));
	y=floor(dist/p)+4;
endfunction
