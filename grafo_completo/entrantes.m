# Uso:
#
#	y=entrantes(i,A)
#
# Donde i debe ser entero y A debe ser una matriz de indices
# y valores de una matriz de incidencia de un grafo.
# Dado i indice de un nodo en un grafo dirigido dado por la
# matriz A, devuelve un vector columna y, con la lista de los
# nodos de A que llegan a i, respetando el sentido en el que
# deben ser recorridas las aristas.
#
function y=entrantes(j,A)
	A=sparse(A(:,1),A(:,2),A(:,3));
	j_1=find((A(:,j)+1)==0); % aristas que terminan en j
	y=find((sum(A(j_1,:),1)-1)==0)'; % nodos contiguos a j
endfunction
