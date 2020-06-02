% -------------------------------------------------------------------
% By JRR May 29, 2020.
% It creates a stochastic matrix by rows
%
% -------------------------------------------------------------------
function P = make_stochastic_matrix(n)
  P = rand(n).*validRandi(n);
  for i=1:n
    P(i,:) = P(i,:)/sum(P(i,:));
  end
end



function R = validRandi(n)
  R = randi([0,1],n);
  while 1    
    if prod(sum(R,1))*prod(sum(R,2))
      return
    else
      R = randi([0,1],n);
    end
  end
end
