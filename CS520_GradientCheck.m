% Project 1: optimization on a sophere
% Gradient Check
clear

global lambda rho N P

lambda = 10 * ones(1,N);
rho = 100;

N = 2;
P = 3;
mu = 0.1;

X = 2*rand(P,N)-1;
X = normalize(X)


% Gradient check for exterior penalty method

% Central difference
p = rand(P, N);
p = 0.001 * normalize(p);
Diff = (obj_ep(X+p,mu)-obj_ep(X,mu))/0.001

% Gradient computation
Gra = zeros(P,N);
for i = 1:N
    for j = 1:N
        if j ~= i
            Gra(:,i) = Gra(:,i) - 2 * (X(:,i)-X(:,j))/(norm(X(:,i)-X(:,j))^4) ;
        end
    end
    Gra(:,i) = Gra(:,i) + 2 * mu * (norm(X(:,i))^2-1) * X(:,i);
end

Gra

Gra = Gra(:)
p = p(:)

transpose(Gra) * p*1000



% Gradient check for augmented Lagrangian system

% Central difference
p = rand(P, N);
p = 0.0001 * normalize(p);
Diff = (al(X+p)-al(X))/0.0001

% Gradient computation
Gra = zeros(P,N);
for i = 1:N
    for j = 1:N
        if j ~= i
            Gra(:,i) = Gra(:,i) - 2 * (X(:,i)-X(:,j))/(norm(X(:,i)-X(:,j))^4);
        end
    end
    Gra(:,i) = Gra(:,i) + 2 * rho * (norm(X(:,i))^2-1) * X(:,i);
    Gra(:,i) = Gra(:,i) + 2 * lambda(i) * X(:,i);
end


Gra = Gra(:);
p = p(:);

transpose(p) * Gra * 10000


function y = al(X)
    global N lambda rho
    y = 0;
    for i = 1:N
        for j = 1:(i-1)
            y = y + 1./(norm(X(:,i)-X(:,j))^2);
        end
        y = y + lambda(i) * (norm(X(:,i))^2-1) + 0.5 * rho * (norm(X(:,i))^2-1)^2;
    end
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



