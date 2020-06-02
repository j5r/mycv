% -------------------------------------------------------------------
% By JRR May 29, 2020.
% It plots the probability distribution of a homogeneous Markov chain 
% since it is given the transition probability matrix P an
% initial distribution pi0 and a number of discrete periods to sample.
%
% -------------------------------------------------------------------
function plotdistrib(P,pi0,n)
  n=n+1;
  pi = zeros(numel(pi0),n);
  pi(:,1) = pi0(:);
  for i=2:n
    pi(:,i) = P' * pi(:,i-1);
  end
  
  plot([0:n-1],pi(1,:),'*-','color',rand(1,3),'lineWidth',2);
  lgs = {'[pi_k]_1'};
  hold on;
  for i=2:size(pi,1)
    plot([0:n-1],pi(i,:),'*-','color',rand(1,3),'lineWidth',2);
    lgs{i}  = ['[pi_k]_',num2str(i)];
  end
  hold off;
  grid on;
  legend(lgs);
  xlabel('tempo k');
  ylabel('Probabilidade de x_k = i');
  %title('Evolução da distribuição pi_0 com matriz P');  
end
