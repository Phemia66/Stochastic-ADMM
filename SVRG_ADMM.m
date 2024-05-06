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
function [xr,E_train,E_test,TM]=SVRG_ADMM(x_0,sample_train,label_train, Num_train,sample_test, label_test,Num_test,...
    A,params)

x_sum = zeros(size(x_0));

if(~exist('params','var'))
    params=struct();
end
    
params = Set_svrg_Params(params);
Num_trail = params.Num_trail;

for T = 1:Num_trail
    %==============================================================================================

    maxoutit = params.outiter;
    max_innerit = params.innerit;
    tol = params.tol;
    eta = params.eta;
    beta = params.beta;
    mu1 = params.mu1;
    mu2 = params.mu2;
    bsize = params.bsize;

    %============================================================
    x_s = x_0;
    count = 2;
    E_train(1,T) = Compute_loss(sample_train,x_s,label_train,Num_train,A,mu1,mu2);
    E_test(1,T) = Compute_loss(sample_test,x_s,label_test,Num_test,A,mu1,mu2);
    
    TM(1,T) = 0;
    
    tic;
 for s = 1:maxoutit
    x = x_s;
    u = A*x;
    z = u;
    
    temp1 = label_train.*(sample_train*x);
    temp2 = -label_train.* exp(temp1)./(1+exp(temp1)).^2;
    v = sample_train'* temp2/Num_train + mu2*x;

    
    %============================================================   
    
      
    md = floor(Num_train/bsize); %number of batch
    
    %============================================================    

    nb=1;
    RD = zeros(1,md*nb);
   

        for j = 1:nb
            idx = randperm(md);
            RD(1,(j - 1)*md + 1:j*md) = idx;
        end

    
    temp_matrix = inv(beta* A'*A + (1/eta)*eye(length(x)));
    
    for i = 1 : max_innerit
        


        
       %%%% update z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Par1 = A*x - u;
        z = Prox(Par1,mu1/beta);

        
       %%% compute stochastic gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            if mod(i, md) == 0
                bdx = md;
            else
                bdx = mod(i,md);
            end

            %%%%% compute sto grdient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            sto = RD(1,bdx);
            Idxsub = (sto - 1)*bsize + 1:sto*bsize;
            temp3 = label_train(Idxsub,:).*(sample_train(Idxsub,:)*x);
            temp4 = -label_train(Idxsub,:) .* exp(temp3)./(1+exp(temp3)).^2;
            sgd1 = sample_train(Idxsub,:)'* temp4 + bsize*mu2*x;
            
            temp5 = label_train(Idxsub,:).*(sample_train(Idxsub,:)*x_s);
            temp6 = -label_train(Idxsub,:) .* exp(temp5)./(1+exp(temp5)).^2;
            sgd2 = sample_train(Idxsub,:)'* temp6 + bsize*mu2*x_s;
            
            Sgd = (sgd1 - sgd2)/bsize + v;
            
            %%%%%%%%%%%%%% update x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
             x = temp_matrix * (beta.* A'*( z + u)+ x./eta - Sgd);
            
            %%%%%%%%%%%%%% update u %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
            u = u - (A*x - z);
            
    end
    
    x_s = x;
    
    tim1=toc;
    errs = Compute_loss(sample_train,x_s,label_train,Num_train,A,mu1,mu2);
    E_train(count,T) = errs;
    E_test(count,T) = Compute_loss(sample_test,x_s,label_test,Num_test,A,mu1,mu2);
    
    TM(count,T) = tim1;

       
       if tim1 >= 10
           break;
       end
        count = count + 1;
        

 end
 

    x_sum = x_sum + x;
   
end

xr = x_sum/Num_trail;
end
        
       

%================================================
% Prox operator.
%================================================
function s = Prox(s,c)
s = sign(s).*(max(abs(s) - c,0));
end
%================================================









