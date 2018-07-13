% Michael R. Buche
% MAE 5730 - Problem 48
% Final Computation Project
% Right-hand side function

function dz = p48_RHS(~,z,q)
    M = vpa(subs(q.M,q.all_vars,z));
    b = vpa(subs(q.b,q.all_vars,z));
    x = M\b;
    dz = zeros(q.num_vars,1);
    for i = 1:q.num_vars
        dz(i) = x(i+q.offset);
    end
end