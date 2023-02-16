%% calculate leaky modes in an approximation sense

thk=[5 5];
dns=[2 2 2];
vs=[500 400 300]; %% for model1
%vs=[400 500 300]; %% for model2

freq=5:80;
vp=2*vs;

cr = modal_v(freq,thk,dns,vs,vp);


%% calculate leaky modes by muller method, theoretical values

[k,cr_t,cr_imag] = leaky_fast(freq,thk,dns,vs,vp);

%% compare them

plot_kai;
