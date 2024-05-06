function s = Compute_time(matrix)
[n1, n2] = size(matrix);
s = zeros(n1,n2);
s(1,:) = matrix(1,:);
for i = 2:n1
    s(i,:) = s(i-1,:) + matrix(i,:);
end
end 