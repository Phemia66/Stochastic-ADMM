function s = Compute_loss(matrix,x,label_vector,Num,A,mu1,mu2)
s = 0;
for i = 1:Num
    s = s + 1/(1 + exp(label_vector(i,1)*matrix(i,:)*x));
end
s = s/Num + mu2*norm(x,2)^2/2 + mu1*norm(A*x,1);
end 