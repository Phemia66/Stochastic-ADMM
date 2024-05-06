function s = SetParams(s, p)
% s = SetDefaultParams(s);
% Sets default parameters
% s: user-specified parameters that are used instead of defaults

if (isfield(s, 'tol') == 0),
    s.tol = 1e-2;
end

if (isfield(s, 'niter') == 0),
    s.niter = 50000; 
end

if (isfield(s, 'epoch') == 0),
    s.epoch=1;
end

if (isfield(s, 'bsize') == 0),
    s.bsize=1000;
end

if (isfield(s,'sigma')==0),
    s.sigma=1;
end

if (isfield(s, 'mu1') == 0),
    s.mu1=1e-5;
end

if (isfield(s, 'mu2') == 0),
    s.mu2=1e-5;
end

if (isfield(s, 'tau') == 0),
    s.tau = 10;
end

if (isfield(s, 'beta') == 0),
    s.beta =1;
end

if (isfield(s, 'eta') == 0),
    s.eta = 0.5;   % for stoADMM
end




    
    
    
    
    
    
    
    
