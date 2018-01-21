# Uso:
#
#	y=a_star(o,d,A)
#
# Donde o y d deben ser enteros y A una matriz de indices de 
# posiciones y valores de una matriz de incidencia de un 
# grafo.
# Dados o y d indices de nodos del grafo descripto por la 
# matriz A, encuentra el camino %*más*) corto entre o y d, en 
# termino de los costos de las aristas, usando un algortimo
# de busqueda de tipo 'A estrella'. Usa como %*función*) 
# %*heurística*) la distancia dada por la norma 1. Considera que
# %*encontró*) una ruta si se encuentra menos de menos de 300 m
# del origen o del destino. Esto debe ser especificado con 
# los signos de o y d.
# Si ambos son positivos, termina si encuentra un camino que
# llegue a menos de 300 m del destino. Si ambos son 
# negativos, termina si encuentra un camino que parta de 
# %*algún*) lugar a menos de 300 m del origen.
# Devuelve un vector fila 'y' con la %*sucesión*) de nodos que 
# tiene el camino que haya encontrado.
# Si origen y destino %*están*) a menos de 300 m, devuelve una
# arista cualquiera
# Si no encuentra %*ningún*) camino, devuelve un vector %*vacío*).
#
# NOTA: usa salientes.m, entrantes.m y arista.m
#
function y=a_star(o,d,A)
	nodos=load('data/input/nodos_amp');
	c_a=load('data/input/costo_aristas');
	closed=[];
	if(o>0) % si el grafo debe recorrerse 'mano'
		if (sum(abs(nodos(o,:)-nodos(d,:)).*[94870 111180])<300) 
		% si origen y destino %*están*) a menos de 300 m
			i_1=(salientes(o,A))(1);
			y=[o i_1];
		else
			h_aux=sum(abs(nodos(o,:)-nodos(d,:)));
			open=[o 0 h_aux];
			came_from=zeros(rows(nodos),1);
			while (rows(open)!=0)
				[x i]=min(open(:,2)+open(:,3));
% guardo el nodo mas barato en current
				current=open(i,:);	
				if (sum(abs(nodos(current(1),:)-nodos(d,:)).*[94870 111180])<300) 
% si estoy a menos de 300 m del destino, termino
					break;
				endif
% agrego el nodo actual a closed
				closed=[closed;open(i,:)];
				i_aux=setdiff([1:rows(open)],i);
% saco el nodo actual de open
				open=open(i_aux,:); 
% busco los vecinos de current
				i_1=salientes(current(1),A); 
				for i=1:rows(i_1)
% si ya esta en closed, sigo
					if (length(find(closed(:,1)==i_1(i))))
						continue;
					endif
					g_aux=current(2)+c_a(find(arista(current(1),i_1(i),A)));
					h_aux=sum(abs(nodos(i_1(i),:)-nodos(d,:)));
					if (length(find(open(:,1)==i_1(i)))) 
% si ya %*está*) en open, lo dejo con el menor costo
						i_aux=find(open(:,1)==i_1(i));
						if(open(i_aux,2)>g_aux)
							open(i_aux,2)=g_aux;
							came_from(open(i_aux,1))=current(1);
						endif
					else	
% si no %*está*) en open, lo agrego
						open=[open;[i_1(i) g_aux h_aux]];
						came_from(i_1(i))=current(1);
					endif
				endfor
			endwhile
			if (rows(open)==0)
				y=[];
			else
% reconstruyo el camino
				flag=current(1);
				y=current(1);
				while (flag!=o)
					y=[came_from(flag);y];
					flag=came_from(flag);
				endwhile
				y=y';
			endif
		endif
% si el grafo debe recorrerse 'contramano'
	else 
		o=-o;
		d=-d;
		if(sum(abs(nodos(o,:)-nodos(d,:)).*[94870 111180])<300) 
			i_1=(entrantes(d,A))(1);
			y=[i_1 d];
		else
			h_aux=sum(abs(nodos(o,:)-nodos(d,:)));
			open=[d 0 h_aux];
			came_from=zeros(rows(nodos),1);
			while (rows(open)!=0)
				[x i]=min(open(:,2)+open(:,3));
				current=open(i,:);
				if (sum(abs(nodos(current(1),:)-nodos(o,:)).*[94870 111180])<300) 
					break;
				endif
				closed=[closed;open(i,:)];
				i_aux=setdiff([1:rows(open)],i);
				open=open(i_aux,:);
				i_1=entrantes(current(1),A);
				for i=1:rows(i_1)
					if (length(find(closed(:,1)==i_1(i))))
						continue;
					endif
					g_aux=current(2)+c_a(find(arista(i_1(i),current(1),A)));
					h_aux=sum(abs(nodos(i_1(i),:)-nodos(o,:)));
					if (length(find(open(:,1)==i_1(i)))) 
						i_aux=find(open(:,1)==i_1(i));
						if(open(i_aux,2)>g_aux)
							open(i_aux,2)=g_aux;
							came_from(open(i_aux,1))=current(1);
						endif	
					else
						open=[open;[i_1(i) g_aux h_aux]];
						came_from(i_1(i))=current(1);
					endif
				endfor
			endwhile
			if (rows(open)==0)
				y=[];
			else
				flag=current(1);
				y=current(1);
				while (flag!=d)
					y=[y;came_from(flag)];
					flag=came_from(flag);
				endwhile
				y=y';
			endif
		endif
	endif
endfunction
