% Project 1: optimization on a sophere
% Solution with augmented Lagrangian system using quasi-Newton and
% Nelder-Mead method
clear
global lambda rho N P

N = 10;
P = 3;
mu = 0.1;

obj_final = 1000;
lambda = 10 * ones(1,N);
rho = 100;
lam_mat = lambda;
for k = 1:1
    X0 = 2*rand(P,N)-1;
    X0 = normalize(X0);
    tic
    % Quasi-Newton method
%     options = optimoptions('fminunc', 'Algorithm', 'quasi-newton','MaxFunctionEvaluations', 20000);
%     [xopt,fval,exitflag,output] = fminunc(@fun2, X0, options);
    % Nelder-Mead simplex algorithm
    options = optimset('MaxFunEvals', 20000, 'MaxIter', 20000);
    [xopt,fval,exitflag,output] = fminsearch(@fun2, X0, options);
    output
    toc 
%     fprintf('xopt: %f, fval: %f\n', xopt, fval);
    disp(['CPU time: ',num2str(toc)]);
    for i = 1:N
        lambda(i) = lambda(i) + rho * (norm(xopt(:,i))^2-1);
    end
    lam_mat = [lam_mat; lambda];
    obj = objective(xopt);
    if obj < obj_final
        obj_final = obj;
        xopt_final = xopt;
    end
end

fprintf('The best answer is: %f',obj_final)

figure(1);
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

figure(2);
plot(lam_mat(:,1), 'linewidth', 2); 
set(gca,'linewidth', 0.75, 'fontsize', 15, 'fontname', 'Arial');
ylabel('value of \lambda_1');
xlabel('iteration');

figure(3);
plot(lam_mat(:,2), 'linewidth', 2); 
set(gca,'linewidth', 0.75, 'fontsize', 15, 'fontname', 'Arial');
ylabel('value of \lambda_2');
xlabel('iteration');

function x = normalize(X)
    [p n] = size(X);
    for i = 1:n
        x(:,i) = X(:,i)/norm(X(:,i));
    end
end

function y = objective(X)
    [p n] = size(X);
    y = 0;
    for i = 1:n
        for j = 1:(i-1)
            y = y + 1/(norm(X(:,i)-X(:,j))^2);
        end
    end
end

function y = fun2(X)
    global N lambda rho
    y = 0;
    for i = 1:N
        for j = 1:(i-1)
            y = y + 1./(norm(X(:,i)-X(:,j))^2);
        end
        y = y + lambda(i) * (norm(X(:,i))^2-1) + 0.5 * rho * (norm(X(:,i))^2-1)^2;
    end
end
