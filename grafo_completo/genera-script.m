A=load(’data/input/matriz_A_2’);
c_a=load(’data/input/costo_aristas’);
c_a=c_a(find(c_a));
p=mean(c_a);
genera_D_G(’data/input/g_OD’,A,p,3,30);
