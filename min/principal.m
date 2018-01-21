global reg_A c_a Delta Gamma
reg_A=load(’data/input/reg_A’);
c_a=load(’data/input/costo_aristas’);
Delta=load(’data/input/Delta’);
Delta=logical(sparse(Delta(:,1),Delta(:,2),Delta(:,3)));
Gamma=load(’data/input/Gamma’);
Gamma=logical(sparse(Gamma(:,1),Gamma(:,2),Gamma(:,3)));
for i=1:rows(Gamma)
	h0(min(find(Gamma(i,:))))=1;
endfor
if (length(h0)<columns(Gamma))
	h0(columns(Gamma))=0;
endif
h0=h0(:);
[h,f,err,k,salida]=min_gradPr(@T,h0,@(x)qp_s(Gamma,x),[100 1e-12],@reg_armijoPr,[100 1000 1 1 0.5 1e -4]);
