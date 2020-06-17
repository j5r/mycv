% RIBEIRO, J. R. June, 12, 2020. jrodrib@usp.br
%
% Filtering of linear dynamic systems [Kumar&Varaiya:Theorem 2.21]
%
% x(k+1) = A * x(k) + G * w(k)                         % state
%   y(k) = C * x(k) + H * v(k)                         % output
% q(k+1) = A * q(k) + L(k+1) * [y(k+1) - C * A * q(k)] % estimator
%
% w and v are independent noise processes with zero mean, 
% and covariance W and V.
%
% x(0) has mean 0 and covariance x_0_covariance
%
%
%%%% #1
%%%% CREATING PARAMETERS
%%%% 
close all,clear,clc;
% dimentions
k_max = 250;
dim_x = 2;
dim_w = 1;
dim_y = 1;  %must to be <= dim_x
dim_v = 1;

% matrices
A = [-5.0027e-4, -2.2887e-2;
      1.6495e-2,  7.6453e-1];

G = 0.1*[0.7906;
         7.6009];
C = [0,1];

H = 1000;

Q = eye(dim_w, dim_w);               % covariance of w(k)
R = eye(dim_v, dim_v);               % covariance of v(k)
x_0_covariance = eye(dim_x, dim_x);  % mean of x(0)


%%%% #2
%%%% ALLOCATING MEMORY
%%%%
% state vector
XXX = zeros(dim_x, k_max);
x_real = XXX;                        % real value of x(k)
x_posteriori = XXX;                  % post-estimative x(k|k)
x_error_posteriori = XXX;            % post-error x(k) - x(k|k)
x_priori = XXX;                      % pre-estimative x(k|k-1)
x_error_priori = XXX;                % pre-error x(k) - x(k|k-1)
clear XXX;

% covariance
SIGMA = zeros(dim_x, dim_x, k_max); 
Sigma_priori = SIGMA;                % pre-covariance Sigma(k|k-1)
Sigma_posteriori = SIGMA;            % post-covariance Sigma(k|k)
L = zeros(dim_x, dim_y, k_max); 
clear SIGMA;

% output y(k)
y = zeros(dim_y, k_max);

% allocating memory for norms
NORM = zeros(1,k_max);
norm_x_real = NORM;
norm_x_posteriori = NORM;
norm_x_priori = NORM;
norm_x_error_posteriori = NORM;
norm_x_error_priori = NORM;
clear NORM;


%%%% #3
%%%% INITIALIZING PARAMETERS
%%%%
% defining an initial value for x(k=1)
x_real(:,1) = randn(dim_x, 1);

% defining an initial value for y(k=1)
y(:,1) = C*x_real(:,1) + H*randn(dim_v,1);

% defining an initial value for L(k=1)
% L  =  P Q^{-1}
% L  = [Q^{-1}' P']'
L(:,:,1) = [    [C*x_0_covariance*C' + H*R*H']'  \...
                [x_0_covariance*C']'                 ]';

% defining an initial value for q(k=1)
x_posteriori(:,1) = L(:,:,1)*y(:,1);

% defining an initial value for estimation error \tilde{x}_{k=1} = x-q
x_error_posteriori(:,1) = x_real(:,1) - x_posteriori(:,1);

% defining an initial value for Sigma_posteriori(1|1)
Sigma_posteriori(:,:,1) = x_0_covariance - L(:,:,1)*C*x_0_covariance;



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
Q = W/tt;
R = V/tt;
clear W V v w r tt;
% Sigma(0) = covariance of x(0)




%%%% #3
%%%% BUILDING THE SYSTEM DYNAMIC
%%%%

for k = 2:k_max
  x_real(:, k) = A*x_real(:, k-1) + G*randn(dim_w, 1);
    
  y(:, k) = C*x_real(:, k) + H*randn(dim_v, 1);  
  
  x_priori(:,k) = A*x_posteriori(:,k-1);
  
  x_error_priori(:,k) = x_real(:,k) - x_priori(:,k);
  
  Sigma_priori(:,:,k) = A*Sigma_posteriori(:,:,k-1)*A' + G*Q*G';
  
  L(:,:,k) = [[C*Sigma_priori(:,:,k)*C' + H*R*H']'\[Sigma_priori(:,:,k)*C']']';
                 
  Sigma_posteriori(:,:,k) = Sigma_priori(:,:,k) - L(:,:,k)*C*Sigma_priori(:,:,k);
  
  x_posteriori(:,k) = x_priori(:,k) + L(:,:,k)*[y(:,k) - C*A*x_posteriori(:,k-1)];
  
  x_error_posteriori(:,k) = x_real(:,k) - x_posteriori(:,k);
  
  
  norm_x_real(k) = norm(x_real(:,k));
  norm_x_priori(k) = norm(x_priori(:,k));
  norm_x_posteriori(k) = norm(x_posteriori(:,k));
  norm_error_priori(k) = norm(x_error_priori(:,k));
  norm_error_posteriori(k) = norm(x_error_posteriori(:,k));
end




figure;
plot(norm_x_real,'-','LineWidth',3,'Color',[0,0,0]);             %black
hold on;
plot(norm_x_priori,'-.','linewidth',2,'Color',[0,0,1]);           %blue
plot(norm_x_posteriori,'-.','linewidth',2,'Color',[1,0,0]);      %red
plot(norm_x_error_priori,'--','linewidth',1.5,'Color',[0,.8,.8]);   %cyan
plot(norm_x_error_posteriori,'-.','linewidth',1,'Color',[.9,0,.9]); %purple
grid on; hold off;
title('Norms ||x||:  x, x(k|k-1), x(k|k), error(k|k-1) and error(k|k)');
legend('x','x(k|k-1)','x(k|k)','error(k|k-1)','error(k|k)');
saveas(1,'norms.pdf');
close all


for pic=1:dim_x
pic2str = num2str(pic);
myfig = figure;
orient(myfig,'rotated');
plot(x_real(pic,:),'-','LineWidth',3,'Color',[0,0,0]);             %black
hold on;
plot(x_priori(pic,:),'-.','linewidth',2,'Color',[0,0,1]);           %blue
plot(x_posteriori(pic,:),'-.','linewidth',2,'Color',[1,0,0]);      %red
plot(x_error_priori(pic,:),'--','linewidth',1.5,'Color',[0,.8,.8]);   %cyan
plot(x_error_posteriori(pic,:),'-.','linewidth',1,'Color',[.9,0,.9]); %purple
grid on; hold off;
title(['[x]_',pic2str,': x, x(k|k-1), x(k|k), error(k|k-1) and error(k|k)']);
legend('x','x(k|k-1)','x(k|k)','error(k|k-1)','error(k|k)');
saveas(pic,[pic2str,'th-coordinate.pdf']);
end


msg= '     mean       std';
disp('G = ')
disp(G)
disp('H = ')
disp(H)
disp('ERROR PRIORI:       x(k) - x(k|k-1)')
disp(msg)
disp([mean(x_error_priori(:,2:end)')',std(x_error_priori(:,2:end)')'])


disp('ERROR POSTERIORI:   x(k) - x(k|k)')
disp(msg)
disp([mean(x_error_posteriori(:,2:end)')',std(x_error_posteriori(:,2:end)')'])

