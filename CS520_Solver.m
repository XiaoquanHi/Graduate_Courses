% Project 1: optimization on a sophere
% Exploring solvers for nonlinear constrianed problems
clear

N = 100;
P = 3;
mu = 0.1;

rng(2)

obj_final = 10000
for k = 1:1
    X0 = 2*rand(P,N)-1;
    X0 = normalize(X0);
    tic
%     option = optimoptions(@fmincon,'Algorithm','interior-point','ConstraintTolerance',1,'OptimalityTolerance',1e-8,'StepTolerance',1e-6,'FunctionTolerance',1e-6,'Display','final-detailed');
%     option = optimoptions(@fmincon,'Algorithm','sqp-legacy','ObjectiveLimit',1e-3,'ScaleProblem','obj-and-constr', 'MaxFunctionEvaluations', 10000);
    option = optimoptions(@fmincon,'Algorithm','active-set','FunctionTolerance',1e-6,'MaxSQPIter',100,'TolConSQP',1e-6, 'MaxFunctionEvaluations', 10000);
    [xopt,fval,exitflag,output,lambda,grad,hessian] = fmincon('fun1',X0,[],[],[],[],[],[],'con1',option);
%     xopt = normalize(xopt);
%     output
%     norm(grad)
    toc 
    disp(['CPU time: ',num2str(toc)]);
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
% x = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
% y = [4.531, 6.691, 1.842, 7.967, 6.024, 3.448, 3.677, 3.201, 2.576];
x = [1e-2, 1e-1, 1];
y = [36.7 6.23 0.25];
semilogx(x,y,'+-','linewidth',2)
set(gca,'linewidth', 0.75, 'fontsize', 15, 'fontname', 'Arial');
xlabel('maximum constraint violation')
ylabel('optimal value of objective function')


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
