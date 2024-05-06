%==============================================================================================
% Stochastic Variance Reduced Gradient ADMM. (SAGA-ADMM)
%==============================================================================================

%==============================================================================================
% Input:                                             | Output:
% x0: Initial image.                                 | xr: Averaged reconstruct image.
% y: Project vector.                                 | F_SVRG_ADMM: record of error.
% paras: struct store the papameter of fan beam.     | TM: record of time.
% Num_trail: number of experiment.
%==============================================================================================
function [xr,E_train,E_test,TM]=stoADMM(x0,sample_train,label_train, Num_train,sample_test, label_test,Num_test,...
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
    tol = params.tol;
    eta = 5;
    beta = 1;
    mu1 = params.mu1;
    mu2 = params.mu2;
    bsize = 500;
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
    md = floor(Num_train/bsize); %number of batch  

    nb=1;
    RD = zeros(1,md*nb);
   
 
        for j = 1:nb
            idx = randperm(md);
            RD(1,(j - 1)*md + 1:j*md) = idx;
        end
     
    %============================================================
    for i = 1 : maxit

        
        
        %%%% update z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Par1 = A*x - u;
        z = Prox(Par1,mu1/beta);

        
        %%% compute stochastic gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     
            if mod(i, md) == 0
                bdx = md;
            else
                bdx = mod(i,md);
            end
            
            sto = RD(1,bdx);
            Idx = (sto - 1)*bsize + 1:sto*bsize;
           
            %%%%% compute sto grdient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp1 = label_train(Idx,1).*(sample_train(Idx,:)*x);
            temp2 = -label_train(Idx,:) .* exp(temp1)./(1+exp(temp1)).^2;
            sgd = sample_train(Idx,:)'* temp2 + bsize*mu2*x;
           
              
     %%%%%%%%%%%%%% update x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
      x = temp_matrix * (beta.* A'*( z + u)+ x./eta - sgd/bsize);
     
     %%%%%%%%%%%%%% update u %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     u = u - (A*x - z);
      
     %========================================================================================================
      % Update of error and time.
     %========================================================================================================
       
    tim1=toc;
    errs = Compute_loss(sample_train,x,label_train,Num_train,A,mu1,mu2);
    E_train(count,T) = errs;
    E_test(count,T) = Compute_loss(sample_test,x,label_test,Num_test,A,mu1,mu2);
    
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
