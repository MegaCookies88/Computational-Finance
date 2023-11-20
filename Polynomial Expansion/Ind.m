function [N] = Ind(m,n)
% Returns position of basis vector
N = (m+n)*(m+n+1)/2 + 1 + n;
end

