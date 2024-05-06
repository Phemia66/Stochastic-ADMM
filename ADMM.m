%==============================================================================================
% Stochastic Variance Reduced Gradient ADMM. (SAGA-ADMM)
%==============================================================================================
function [xr,E_train,E_test,TM]=ADMM(x0,sample_train,label_train, Num_train,sample_test, label_test,Num_test,...
    A,params)

x_sum = zeros(size(x0));

if(~exist('params','var'))
    params=struct();
end
    
params = SetParams(params);
epoch = params.epoch;


for T = 1:epoch
    %==============================================================================================

    maxit = params.niter;
    tol = 1e-2;
    eta = 0.5;
    beta = 3;
    mu1 = params.mu1;
    mu2 = params.mu2;

    
    
    %============================================================
   
    u = A*x0;
    y = u;
    x = x0;
    
    E_train(1,T) = Compute_loss(sample_train,x,label_train,Num_train,A,mu1,mu2);
    E_test(1,T) = Compute_loss(sample_test,x,label_test,Num_test,A,mu1,mu2);
    TM(1,T) = 0;
    %============================================================    

    temp_matrix = inv(beta* A'*A + (1/eta)*eye(length(x)));
    
    count = 2;
    tic;
    %============================================================
    for i = 1 : maxit
  
        
        
        %%%% update z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Par1 = A*x - u;

        z = Prox(Par1,mu1/beta);

        
        %%% compute stochastic gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           
            %%%%% compute sto grdient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp1 = label_train.*(sample_train*x);
            temp2 = -label_train.* exp(temp1)./(1+exp(temp1)).^2;
            gd_full = sample_train'* temp2/Num_train + mu2*x;
          
           
              
     %%%%%%%%%%%%%% update x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
      x = temp_matrix * (beta.* A'*( z + u)+ x./eta - gd_full);
     
     %%%%%%%%%%%%%% update u %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     u = u - (A*x - z);
      
     %========================================================================================================
      % Update of error and time.
     %========================================================================================================
       
        
    tim1=toc;
    errs = Compute_loss(sample_train,x,label_train,Num_train,A,mu1,mu2);
    E_train(count,T) = errs;
    E_test(count,T) = Compute_loss(sample_test,x,label_test,Num_test,A,mu1,mu2);
    E_acc(count, T) = Predict_test(sample_test, label_test, x, Num_test);
    
    TM(count,T) = tim1;
 
       
       if tim1 >= 10
           break;
       end
        
        count = count + 1;
        

        %========================================================================================================
    end
    

    x_sum = x_sum + x;

end

xr = x_sum/epoch;
end

 

%================================================
% Prox operator.
%================================================
function s = Prox(s,c)
s = sign(s).*(max(abs(s) - c,0));
end
%================================================







