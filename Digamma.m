function d = Digamma(x,n)

v = zeros(1,n);

for m = 1:n
    v(m) = 1/m - 1/(m+x-1);
end

d = psi(1) + sum(v);

end