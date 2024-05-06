function s = Set_sarah_Params(s, p)


if (isfield(s, 'tol') == 0),
    s.tol = 1e-2;
end

if (isfield(s, 'Num_trail') == 0),
    s.Num_trail = 1;
end

if (isfield(s, 'outiter') == 0),
    s.outiter = 30000; 
end

if (isfield(s, 'innerit') == 0),
    s.innerit = 10; 
end

if (isfield(s, 'epoch') == 0),
    s.epoch=1;
end

if (isfield(s, 'bsize') == 0),
    s.bsize=100;
end

if (isfield(s,'sigma')==0),
    s.sigma=1;
end

if (isfield(s, 'tau') == 0),
    s.tau = 0.5;
end

if (isfield(s, 'beta') == 0),
    s.beta =1;
end

if (isfield(s, 'mu1') == 0),
    s.mu1=1e-5;
end

if (isfield(s, 'mu2') == 0),
    s.mu2=1e-5;
end



