function [c,ceq] = con1(X)
    c = [];
    ceq = X'*X-1;
end