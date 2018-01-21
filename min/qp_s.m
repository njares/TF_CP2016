# Uso:
#
#	x1=qp_s(Gamma,h0)
#
# Donde Gamma debe ser una matriz y h0 un vector
# Toma la matriz Gamma de gxr, que es una matriz en bloques,
# donde todos los bloques tiene 1 fila, los bloques de la
# diagonal tienen todas sus entradas iguales a 1, y el resto
# son nulos, y el vector h0 de rx1 y realiza la proyeccion de
# h0 en el conjunto:
#
#       Omega = {h | Gamma*h=ones(g,1) , h >= 0 }
#
# Resolviendo subproblemas que aprovechan la estrucutra de
# Gamma.
# Devuelve el vector x1 en Omega, que es la proyeccion de h0
#
# NOTA: usa qp, funcion provista por Octave
#
function x1=qp_s(Gamma,h0)
	x1=h0;
	j=1;
	for i=1:rows(Gamma)
		m=length(find(Gamma(i,:)));
		h_aux=h0([j:j+m-1]);
		x0=(1/m)*ones(m,1);
		x1([j:j+m-1])=qp(x0,speye(m),-h_aux,ones(1,m),1,zeros(m,1),[]);
		j=j+m;
	endfor
endfunction
