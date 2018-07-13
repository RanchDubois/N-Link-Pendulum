% Michael R. Buche
% MAE 5730 - Problem 48
% Final Computation Project
% Simple automatic parameters function

function p = p48_simple_system(p)

    % Uniform links of unit mass and unit length
    p.g = 9.81; % do not make zero
    p.l = ones(1,p.N);
    p.d = p.l/2;
    p.m = ones(1,p.N);
    p.I_G = p.m.*p.l.^2/12;
    
    % Links start at rest in line at 45 degrees
    if p.nice_initial_conditions
        theta0 = pi/4;
        p.icv = [theta0*ones(1,p.N) zeros(1,p.N)];
        
    % Links start at rest in an upward zig-zag
    else
        angle_1 = pi/2;
        angle_2 = pi;
        half_1 = angle_1*toeplitz(mod(1,2),mod(1:p.N,2));
        half_2 = angle_2*(1 - toeplitz(mod(1,2),mod(1:p.N,2)));
        p.icv = [half_1+half_2 zeros(1,p.N)];
    end
    
    % Torque applied at each hinge
    if p.torques
        p.T = ones(1,p.N);
    end
    
    % Linear torsional springs at each hinge
    if p.springs
        p.k = ones(1,p.N);
    end
    
    % Linear angular viscous friction at each hinge
    if p.friction
        p.c = ones(1,p.N);
    end
end