function [f,Df,D2f]=T(h,der)
% ejecutar [f,Df,D2f]=T(x,der)
% donde T(h)= costo 
%  der(1)=0 no calcula f(x) (default der(1)=1),
%  der(2)=0 no calcula Df(x) (default der(2)=0),
%  der(3)=0 no calcula D2f(x) (default der(3)=0).
f=[]; Df=[]; D2f=[];
d=[true;false;false];
if nargin>1, d(1:length(der))=der; end
h=h(:);
global reg_A c_a Delta
v=Delta'*h;
if d(1)
	w=v;
	w(1:reg_A(1))=v(1:reg_A(1))-exp(-v(1:reg_A(1)));
	f=sum(c_a.*w);
end
if d(2)
	z=[1+exp(-v(1:reg_A(1)));ones(length(v)-reg_A(1),1)];
	z=c_a.*z;
	Df=Delta*z;
end
if d(3)
  D2f=[];
end

