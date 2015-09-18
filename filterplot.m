
function [D] = filterplot(X)

X=X';

[m n] = size(X);


w = round(sqrt(n));
h = (n / w);

c = floor(sqrt(m));
r = ceil(m / c);

p = 1;

D = - ones(p + r * (h + p),p + c * (w + p));

k = 1;
for j = 1:r
    for i = 1:c
        D(p + (j - 1) * (h + p) + (1:h), p + (i - 1) * (w + p) + (1:w)) = reshape(X(k, :), [h, w]) / max(abs(X(k, :)));
        k = k + 1;
    end
    
end

end