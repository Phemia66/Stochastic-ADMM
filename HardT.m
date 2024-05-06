function c=HardT(alpha,beta,lambda)
% c=alpha;
    
% c=double((alpha.^2-2*lambda/beta)>0).*alpha;
c = alpha .* (alpha > sqrt(2*lambda/beta));

end