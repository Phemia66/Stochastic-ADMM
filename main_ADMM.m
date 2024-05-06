%============================================================================
%===========================================================================


clc;
%===========================================================================
x0= zeros(n + 1,1);
%x0= rand(n + 1,1);
%===========================================================================
%===========================================================================
disp('Computation of SAGA-ADMM')
params.tol =0.01;

[x_ADMM_saga,Etr_ADMM_saga,Ets_ADMM_saga,Tim_ADMM_saga] = SAGA_ADMM(x0,sample_train,label_train, Num_train, sample_test, label_test, Num_test,A,params);

Ets_ADMM_saga(length(Ets_ADMM_saga))

%===========================================================================
disp('Computation of SARAH-ADMM')
params.tol =0.01;

[x_ADMM_sah,Etr_ADMM_sah,Ets_ADMM_sah,Tim_ADMM_sah] = SARAH_ADMM(x0,sample_train,label_train, Num_train, sample_test, label_test, Num_test,A,params);

Ets_ADMM_sah(length(Ets_ADMM_sah))

%===========================================================================
disp('Computation of SVRG-ADMM')
params.tol =0.01;

[x_ADMM_svrg,Etr_ADMM_svrg,Ets_ADMM_svrg,Tim_ADMM_svrg] = SVRG_ADMM(x0,sample_train,label_train, Num_train, sample_test, label_test, Num_test,A,params);

Ets_ADMM_svrg(length(Ets_ADMM_svrg))

%===========================================================================
disp('Computation of Stochastic-ADMM')
params.tol =0.01;

[x_ADMM_sto,Etr_ADMM_sto,Ets_ADMM_sto,Tim_ADMM_sto] = stoADMM(x0,sample_train,label_train, Num_train, sample_test, label_test, Num_test,A,params);

Ets_ADMM_sto(length(Ets_ADMM_sto))

%===========================================================================
disp('Computation of ADMM')
params.tol =0.01;

[x_ADMM,Etr_ADMM,Ets_ADMM,Tim_ADMM] = ADMM(x0,sample_train,label_train, Num_train, sample_test, label_test, Num_test,A,params);

Ets_ADMM(length(Ets_ADMM))



%===========================================================================
figure()
plot(Tim_ADMM_saga(1:30:length(Tim_ADMM_saga)),Ets_ADMM_saga(1:30:length(Ets_ADMM_saga)),'-- m*','LineWidth',1);
hold on;
plot(Tim_ADMM_sah(1:30:length(Tim_ADMM_sah)),Ets_ADMM_sah(1:30:length(Ets_ADMM_sah)),'-- g>','LineWidth',1);
hold on;
plot(Tim_ADMM_svrg(1:30:length(Tim_ADMM_svrg)),Ets_ADMM_svrg(1:30:length(Ets_ADMM_svrg)),'-- b+','LineWidth',1);
hold on;
plot(Tim_ADMM_sto(1:30:length(Tim_ADMM_sto)),Ets_ADMM_sto(1:30:length(Ets_ADMM_sto)),'-- co','LineWidth',1);
hold on;
plot(Tim_ADMM(1:30:length(Tim_ADMM)),Ets_ADMM(1:30:length(Ets_ADMM)),'-- r^','LineWidth',1);
legend('SAGA','SARAH','SVRG','SADMM','ADMM','Location','southoutside','Orientation','horizontal');
legend('boxoff')
title('a8a-test-loss');
xlabel('CPU-time (second)');
ylabel('Test loss');









