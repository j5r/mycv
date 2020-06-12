% RIBEIRO, J. R. June, 12, 2020. jrodrib@usp.br
%
% Filtering of linear dynamic systems with exogenous input
%
% x(k+1) = A * x(k) + G * w(k) + m(k)          % state
%   y(k) = C * x(k) + H * v(k)                 % output
% q(k+1) = A * q(k) + L(k) * [y(k) - C * q(k)] % estimator
%
% w and v are independent noise processes with zero mean, 
% and covariance W and V.
%

%%%% #1
%%%% CREATING PARAMETERS
%%%% 
close all,clear,clc;
% dimentions
k_max = 250;
dim_x = 5;
dim_w = 3;
dim_y = 5;
dim_v = 5;

% matrices
A = eye(dim_x, dim_x);
G = eye(dim_x, dim_w);
C = eye(dim_y, dim_x);
H = eye(dim_y, dim_v);

% exogenous inputs
m = randn(dim_x, k_max);

% mean of x(0)
bar_x_0 = zeros(dim_x, 1);

% for uncomment, do %%{ and for comment, do %{
%{
A = randn(dim_x, dim_x);
A = 0.1*A/max(abs(eig(A)));
G = 0.1*randn(dim_x, dim_w);
C = randn(dim_x, dim_x);
while 1
H = randn(dim_x, dim_v);
  if min(eig(H*H'))>0
    break
  end
end
H = 0.001*H;
%}

%%%% #2
%%%% ALLOCATING MEMORY AND INITIALIZING PARAMETERS
%%%%

% allocating memory and giving an initial value for x(k)
x_real = zeros(dim_x, k_max);
x_real(:,1) = bar_x_0 + randn(dim_x, 1);

% allocating memory and giving an initial value for q(k)
x_estimated = zeros(dim_x, k_max);
x_estimated(:,1) = bar_x_0;

% allocating memory and defining an initial value for estimation errorr
x_error = zeros(dim_x, k_max);
x_error(:,1) = x_real(:,1) - x_estimated(:,1);

% allocating memory and defining an initial value for mean estimation errorr
x_error_mean = zeros(dim_x, k_max);
% x_error_mean(:,1) is naturally null

% output y(k)
y = zeros(dim_y, k_max);
y(:,1) = C*x_real(:, 1) + H*randn(dim_v, 1);


%%%%%%%%% calculating covariance of x(0), W, and V
Sigma = zeros(dim_x, dim_x, k_max);
tt = 1e4;
W = zeros(dim_w);
V = zeros(dim_v);
for i=1:tt
  r = randn(dim_x, 1);
  w = randn(dim_w, 1);
  v = randn(dim_v, 1);
  Sigma(:,:,1) = Sigma(:,:,1) + r*r';
  W = W + w*w';
  V = V + v*v';
end
Sigma(:,:,1) = Sigma(:,:,1)/tt;
W = W/tt;
V = V/tt;
clear v w r tt;
% Sigma(0) = covariance of x(0)

% trace of Sigma
trace_Sigma = zeros(1,k_max);
trace_Sigma(1) = trace(Sigma(:,:,1));

% L(0)
L = zeros(dim_x, dim_y, k_max);
L(:,:,1) = [[C*Sigma(:,:,1)*C'+H*V*H']\[m(:,1)*x_error_mean(:,1)'*C'-...
        A*Sigma(:,:,1)*C']']';


% norms
norm_x_real = zeros(1, k_max);
norm_x_real = norm(x_real(:,1));

norm_x_estimated = zeros(1, k_max);
norm_x_estimated = norm(x_estimated(:,1));

norm_x_error = zeros(1, k_max);
norm_x_error = norm(x_error(:,1));

norm_x_error_mean = zeros(1, k_max);
norm_x_error_mean = norm(x_error_mean(:,1));

%%%% #3
%%%% BUILDING THE SYSTEM DYNAMIC
%%%%

for k = 2:k_max
  % x(k)
  x_real(:, k) = A*x_real(:, k-1) + G*randn(dim_w, 1) + m(:,k-1);
  
  %y(k)
  y(:, k) = C*x_real(:, k) + H*randn(dim_v, 1);  
  
  %q(k)
  x_estimated(:, k) = A*x_estimated(:, k-1) + ...
              L(:,:,k-1)*(y(:, k-1) - C*x_estimated(:, k-1));
  %errorr(k)
  x_error(:,k) = x_real(:, k) - x_estimated(:, k);
  
  %mean_errorr(k)
  x_error_mean(:,k) = (A-L(:,:,k-1)*C) * x_error_mean(:,k-1) + m(:,k-1);
  
  %second_moment_errorr(k)
  Sigma(:,:,k) = (A-L(:,:,k-1)*C)*Sigma(:,:,k-1)*(A - L(:,:,k-1)*C)' +...
        L(:,:,k-1)*H*V*H'*L(:,:,k-1)' + G*W*G' + m(:,k-1)*m(:,k-1)' +...
        (A-L(:,:,k-1)*C)*x_error_mean(:,k-1)*m(:,k-1)' +...
        m(:,k-1)*x_error_mean(:,k-1)'*(A-L(:,:,k-1)*C)';
  
  %L(k)
  L(:,:,k) = [[C*Sigma(:,:,k)*C'+H*H']'\[m(:,k)*x_error_mean(:,k)'*C'+...
        A*Sigma(:,:,k)*C']']';
  
  % norms
  trace_Sigma(k) = trace(Sigma(:,:,k));
  
  norm_x_real(k) = norm(x_real(:,k));
  
  norm_x_estimated(k) = norm(x_estimated(:,k));
  
  norm_x_error(k) = norm(x_error(:,k));
  
  norm_x_error_mean(k) = norm(x_error_mean(:,k));
  
end
%figure;
%plot(trace_Sigma,'b-*')
%title('estimation quadratic error')

figure;
plot(norm_x_real,'b-','linewidth',2); hold on;
plot(norm_x_estimated,'r-.','linewidth',2.2);
plot(norm_x_error,'m-.','linewidth',2);
plot(norm_x_error_mean,'g-.','linewidth',2);
grid on; hold off;
title('Norms: real, estimated, error and expected error')
legend('real','estimated','error','expected error');




