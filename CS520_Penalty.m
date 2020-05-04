% Project 1: optimization on a sophere
% Penalty method with backtracking line search
clear

N = 30;
P = 3;
mu = 0.1;

obj_final = 1000000
for k = 1:1000
    X = 2*rand(P,N)-1;
    X = normalize(X);
    xopt = Ex_pen(X,mu);
    xopt = normalize(xopt);
    obj = objective(xopt);
    if obj < obj_final
        obj_final = obj;
        xopt_final = xopt;
    end
end


fprintf('The best answer is: %f',min(obj_final))

figure(2);
sphere(50);
[x y z]=sphere();
s = surf(1*x,1*y,1*z, 'FaceAlpha',0.5);
axis equal;
s.EdgeColor = 'none';
hold on 

x = xopt_final(1,:);
y = xopt_final(2,:);
z = xopt_final(3,:);
scatter3(x,y,z,'filled')


% Exterior penalty method
function xopt = Ex_pen(X,mu)
    [p n] = size(X);
    iteration = 0;
    eps = 1e-6;
    error = 10;
    err_vec = error;
    xopt = X;
    while error > eps && iteration < 1000
        iteration = iteration + 1;
        x_old = xopt;
        grad = zeros(p,n);
        for i = 1:n
            for j = 1:(i-1)
                grad(:,i) = grad(:,i) - 4*(X(:,i)-X(:,j))/(norm(X(:,i)-X(:,j))*norm(X(:,i)-X(:,j))*norm(X(:,i)-X(:,j))*norm(X(:,i)-X(:,j)));
            end
            grad(:,i) = grad(:,i) + mu*(norm(X(:,i))^2-1)*X(:,i);
        end
        alpha = backtr(0.5,x_old,grad,1e-4,0.5,1e-8,mu);
        xopt = x_old + alpha.*grad;
        error = norm(xopt-x_old);
        fprintf('iteration: %d, error: %d \n', iteration, error)
%         display(xopt);
        err_vec = [err_vec error];
    end
%     figure(1);
%     plot(log(err_vec));
%     ylabel('log(||x_{k+1}-x_k||)');
%     xlabel('iteration');
end

function x = normalize(X)
    [p n] = size(X);
    for i = 1:n
        x(:,i) = X(:,i)/norm(X(:,i));
    end
end

function y = obj_ep(X,mu)
    [p n] = size(X);
    y = 0;
    for i = 1:n
        for j = 1:(i-1)
            y = y + 1./(norm(X(:,i)-X(:,j)).*norm(X(:,i)-X(:,j)));
        end
        y = y + 0.25*mu*(norm(X(:,i)).^2-1).^2;
    end
end

function y = objective(X)
    [p n] = size(X);
    y = 0;
    for i = 1:n
        for j = 1:(i-1)
            y = y + 1./(norm(X(:,i)-X(:,j))^2);
        end
    end
end


% BACKTRACKING ARMIJO-TYPE
function [alpha] = backtr(alpha_guess,Xk,dk,gamma,delta,rhok,mu)
% INPUT:(*) alpha _guess - current steplength (1*1) [>0];
%       (*) Xk           - current iterate    (N*1);
%       (*) dk           - search direction   (N*1);
%           gamma        - constant provided by the user (1*1) [>0];
%           delta        - constant provided by the user (1*1) into the range [0,  1];
%           rhok         - constant provided by the user (1*1) into the range [0,  1];
%       (*) F            - function handle of the objective function (RN->R );
% OUTPUT:    alpha - value of alpha whether the condition holds (1*1);
    % positive direction (+)alpha
    alpha = alpha_guess;
    while (obj_ep(Xk+alpha.*dk,mu)>obj_ep(Xk,mu)-gamma*alpha^2*(norm(dk))^2) 
        if (alpha*norm(dk) < rhok)   
            alpha  = 0;              % <-- failure to search for a value of alpha nonzero
        else
            alpha = alpha*delta;     % <-- reduction of the steplength
        end
    end 
    alpha1 = alpha;
    F1     = obj_ep(Xk+alpha1.*dk, mu)-(obj_ep(Xk, mu)-gamma*alpha1^2*(norm(dk))^2);    
    % negative direction (-)alpha
    alpha = alpha_guess;
    while (obj_ep(Xk-alpha.*dk, mu)>obj_ep(Xk, mu)-gamma*alpha^2*(norm(dk))^2)  
        if (alpha*norm(dk) < rhok)
            alpha   = 0;              % <-- failure to search for a value of alpha nonzero
        else
            alpha = alpha*delta;      % <-- reduction of the steplength
        end
    end
    alpha2 = -alpha;
    F2     = obj_ep(Xk+alpha2.*dk,mu)-(obj_ep(Xk,mu)-gamma*alpha2^2*(norm(dk))^2);
    % choice of the value of alpha for which it is provided with sufficient reduction 
    if (F1<F2)           
        alpha = alpha1;
    else
        alpha = alpha2;
    end  
    
end


