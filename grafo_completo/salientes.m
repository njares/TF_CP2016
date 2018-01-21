# Uso:
#
#	y=salientes(i,A)
#
# Donde i debe ser entero y A debe ser una matriz de indices
# y valores de una matriz de incidencia de un grafo.
# Dado i indice de un nodo en un grafo dirigido dado por la
# matriz A, devuelve un vector columna y, con la lista de los
# nodos de A a los que se puede llegar desde i, respetando el
# sentido en el que deben ser recorridas las aristas.
#
function y=salientes(i,A)
	A=sparse(A(:,1),A(:,2),A(:,3));
	i_1=find((A(:,i)-1)==0); % aristas que se inician en i
	y=find((sum(A(i_1,:),1)+1)==0)'; % nodos contiguos a i
endfunction
