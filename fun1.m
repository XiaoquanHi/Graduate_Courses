function y = fun1(X)
    [p n] = size(X);
    y = 0;
    for i = 1:n
        for j = 1:(i-1)
            y = y + 1./(norm(X(:,i)-X(:,j))^2);
        end
    end
end