function Aspec = Spectraldecomposition(A)

[rV,D,lV] = eig(A);

Aspec=zeros(length(A));
for i=1:rank(A)
    Aspec= Aspec + D(i,i) * ((rV(:,i)) * lV(:,i)');
end

disp(Aspec)

end






