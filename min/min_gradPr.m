function [x,f,err,k,salida]=min_gradPr(fun,x,PrX,parmet,funregla,parregla)
% ejecutar [x,f,err,k,salida]=min_cauchy(fun,x0,PrX)
%        o [x,f,err,k,salida]=min_cauchy(fun,x0,PrX,parmet)
%        o [x,f,err,k,salida]=min_cauchy(fun,x0,PrX,parmte,funregla)
%        o [x,f,err,k,salida]=min_cauchy(fun,x0,PrX,parmte,funregla,parregla)
%
% donde parmet(1)=cantidad maxima de iteraciones (default=100)
%       parmet(2)=tolerancia de error del gradiente (default=1e-7)
%       funregla=funcion regla de busqueda lineal (default=@reg_armijo)
%       parregla=parametros para la regla escogida
%
% salida=0, gradiente menor que tolerancia.
% salida=1, maximo de iteraciones alcanzado.
% salida=3, maximo de iteraciones en busqueda lineal.
pm=[100,1e-7]; regla=@reg_armijoPr; pr=100;
if nargin>3, pm(1:length(parmet))=parmet; end
if nargin>4&&~isempty(funregla), regla=funregla; end
if nargin==6, pr(1:length(parregla))=parregla; end
maxit=pm(1); tol=pm(2);
x=PrX(x);
[f,df]=fun(x,[1,1]);
err=norm(x-PrX(x-df));
for k=0:maxit
	printf ("iter_min=%d fun=%e err=%e\n",k,f,err)
  if err<tol
    salida=0;
    return
  end
  [x,t,ki]=regla(fun,x,PrX,f,df,pr);
  if ki>=pr(1)
    salida=3;
    return
  end
  [f,df]=fun(x,[1,1]);
  err=norm(x-PrX(x-df));
end
salida=1;
