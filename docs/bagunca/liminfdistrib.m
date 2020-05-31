% -------------------------------------------------------------------
% By JRR May 29, 2020.
% Solves the problem, for a stochastic matrix P (by rows)
%
%         pi_inf = P' pi_inf
%
% The output is all stationary distributions
% Warning! These distributions may depend upon the initial state pi_0
%
% -------------------------------------------------------------------
function result = liminfdistrib(P,tolerance)
  if nargin ==  1
    tolerance = eps(10);
  end
  [eigenvecs,eigenvals] = eig(P');
  where_eigenval_is_1 = abs(real(diag(eigenvals))-1) < tolerance; % logical value
  result = eigenvecs( :, where_eigenval_is_1);  
  for i=1:size(result,2)
      result(:,i) = result(:,i)/sum(result(:,i));
  end
end
