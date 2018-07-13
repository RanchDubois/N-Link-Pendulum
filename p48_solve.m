% Michael R. Buche
% MAE 5730 - Problem 48
% Final Computation Project
% Solve for the motion

function [t,fx_G,fy_G,fx_end,fy_end,ftheta,fdtheta,fv_G] = p48_solve(tspan,p)

    % Start timer to record time taken to solve
    tic
    
    % Halt certain cases that are not permissible or not yet set up
    if strcmp(p.method,'minimal_Lagrange')
        assert(~p.friction,'Cannot use Lagrange equations while there is friction!')
        assert(~p.torques,'Lagrange equations method not set up for external torques.')
    end
    if strcmp(p.method,'maximal')
        assert(~p.friction,'Maximal coordinates method not set up for friction.')
        assert(~p.torques,'Maximal coordinates method not set up for external torques.')
        assert(~p.springs,'Maximal coordinates method not set up for torsional springs.')
    end
    
    % Display some information about the system
    disp(['Intializing variables for ',num2str(p.N), ' link pendulum,'])
    if p.simple_system
        if p.nice_initial_conditions
            disp('    nice initial conditions,')
        else
            disp('    not nice initial conditions,')
        end
        if p.torques
            disp('    torques applied at hinges,')
        else
            disp('    no applied torques,')
        end
        if p.springs
            disp('    torsional springs included,')
        else
            disp('    no torsional springs,')
        end
        if p.friction
            disp('    hinge frictions included,')
        else
            disp('    no hinge frictions,')
        end
    end

    % Fixed frame unit vectors
    i_hat = [1 0 0];
    j_hat = [0 1 0];
    k_hat = [0 0 1];

    % Minimal coordinates method
    if strcmp(p.method,'minimal_AMB')
        disp('    using minimal coordinates and AMB...')

        % Initialize symbolic variables
        syms t 
        theta = sym([]);
        er = sym([]);
        et = sym([]);
        rG_O = sym([]);
        rG_A = sym([]);
        rA_B = sym([]);
        r_hingei_Gj = sym([]);
        aG_O = sym([]);
        aG_A = sym([]);
        aA_B = sym([]);
        M_A = sym([]);
        H_dot_A = sym([]);
        eqns = sym([]);

        % Distance and acceleration vectors
        disp('Creating distance and acceleration vectors...')
        for i = 1:p.N
            theta(i) = symfun(sym(sprintf('theta_%d(t)', i)), t);
            er(i,:) = [cos(theta(i)) sin(theta(i)) 0];
            et(i,:) = cross(k_hat,er(i,:));
            rA_B(i,:) = p.l(i)*er(i,:);
            rG_A(i,:) = p.d(i)*er(i,:);
            rG_O(i,:) = rG_A(i,:);
            if (i > 1)
                for j = 1:(i-1)
                    rG_O(i,:) = rG_O(i,:) + rA_B(j,:);
                end
            end
            aA_B(i,:) = cross(diff(theta(i),2)*k_hat,rA_B(i,:))...
                - diff(theta(i))^2*rA_B(i,:);
            aG_A(i,:) = cross(diff(theta(i),2)*k_hat,rG_A(i,:))...
                - diff(theta(i))^2*rG_A(i,:);
            aG_O(i,:) = aG_A(i,:);
            if (i > 1)
                for j = 1:(i-1)
                    aG_O(i,:) = aG_O(i,:) + aA_B(j,:);
                end
            end
        end

        % Distance from ith hinge to jth center of mass
        disp('Creating additional distance vectors...')
        for i = 1:p.N
            for j = 1:p.N
                r_hingei_Gj(i,:,j) = rG_A(j,:);
                if (j > 1)
                    for k = i:(j-1)
                        r_hingei_Gj(i,:,j) = r_hingei_Gj(i,:,j) + rA_B(k,:);
                    end
                end
            end
        end

        % Angular momentum balance equations
        disp('Writing equations of motion...')
        for i = 1:p.N
            M_A(i,:) = cross(rG_A(i,:),p.m(i)*p.g*i_hat);
            H_dot_A(i,:) = p.m(i)*cross(rG_A(i,:),aG_O(i,:))...
                + p.I_G(i)*diff(theta(i),2)*k_hat;
            if (i < p.N)
                for j = (i+1):p.N
                    M_A(i,:) = M_A(i,:) + cross(r_hingei_Gj(i,:,j),p.m(j)*p.g*i_hat);
                    H_dot_A(i,:) = H_dot_A(i,:)...
                        + p.m(j)*cross(r_hingei_Gj(i,:,j),aG_O(j,:))...
                        + p.I_G(j)*diff(theta(j),2)*k_hat;
                end
            end
            if p.torques
                M_A(i,3) = M_A(i,3) + p.T(i);
            end
            if p.springs
                if i == 1
                    M_A(1,3) = M_A(1,3) - p.k(1)*theta(1);
                else
                    M_A(i,3) = M_A(i,3) - p.k(i)*(theta(i) - theta(i-1));
                end
            end
            if p.friction
                if i == 1
                    M_A(1,3) = M_A(1,3) - p.c(1)*diff(theta(1));
                else
                    M_A(i,3) = M_A(i,3) - p.c(i)*(diff(theta(i)) - diff(theta(i-1)));
                end
            end
            eqns(i) = M_A(i,3) - H_dot_A(i,3);
        end
    end

    % Maximal coordinates method
    if strcmp(p.method,'maximal')
        disp('    using maximal coordinates...')

        % Initialize symbolic variables
        syms t 
        theta = sym([]);
        x_G = sym([]);
        y_G = sym([]);
        x_end = sym(zeros(1,p.N));
        y_end = sym(zeros(1,p.N));
        R_x = sym([]);
        R_y = sym([]);
        eqns_LMB = sym([]);
        eqns_AMB = sym([]);
        eqns_con = sym([]);

        % Symbolic functions of time
        for i = 1:p.N
            theta(i) = symfun(sym(sprintf('theta_%d(t)', i)), t);
            x_G(i) = symfun(sym(sprintf('x_G_%d(t)', i)), t);
            y_G(i) = symfun(sym(sprintf('y_G_%d(t)', i)), t);
            R_x(i) = symfun(sym(sprintf('R_x_%d(t)', i)), t);
            R_y(i) = symfun(sym(sprintf('R_y_%d(t)', i)), t);
        end

        % Linear and angular momentum balance and constraint equations
        disp('Writing equations of motion')
        disp('    and constraint equations...')
        for i = 1:p.N
            if (i == 1)
                eqns_LMB(1,1) = p.m(1)*diff(x_G(1),2) - p.m(1)*p.g - R_x(1) - R_x(2);
                eqns_LMB(1,2) = p.m(1)*diff(y_G(1),2) - R_y(1) - R_y(2);
                eqns_AMB(1) = -p.d(1)*(R_x(1)*sin(theta(1)) - R_y(1)*cos(theta(1)))...
                    - (p.l(1) - p.d(1))*(R_y(2)*cos(theta(1)) - R_x(2)*sin(theta(1)))...
                    + p.I_G(1)*diff(theta(1),2);
                x_end(1) = p.l(1)*x_G(1)/p.d(1);
                y_end(1) = p.l(1)*y_G(1)/p.d(1);
                eqns_con(1,1) = diff(x_G(1) - p.d(1)*cos(theta(1)),2);
                eqns_con(1,2) = diff(y_G(1) - p.d(1)*sin(theta(1)),2);
            else
                if (i < p.N)
                    eqns_LMB(i,1) = p.m(i)*diff(x_G(i),2) - p.m(1)*p.g...
                        + R_x(i) - R_x(i+1);
                    eqns_LMB(i,2) = p.m(i)*diff(y_G(i),2) + R_y(i) - R_y(i+1);
                    eqns_AMB(i) = -p.d(i)*(R_y(i)*cos(theta(i))...
                        - R_x(i)*sin(theta(i)))...
                        - (p.l(i) - p.d(i))*(R_y(i+1)*cos(theta(i))...
                        - R_x(i+1)*sin(theta(i)))...
                        + p.I_G(i)*diff(theta(i),2);
                else
                    eqns_LMB(p.N,1) = p.m(p.N)*diff(x_G(p.N),2) - p.m(1)*p.g + R_x(p.N);
                    eqns_LMB(p.N,2) = p.m(p.N)*diff(y_G(p.N),2) + R_y(p.N);
                    eqns_AMB(p.N) = -p.d(p.N)*(R_y(p.N)*cos(theta(p.N))...
                        - R_x(p.N)*sin(theta(p.N)))...
                        + p.I_G(p.N)*diff(theta(p.N),2);
                end
                x_end(i) = x_end(i-1) + p.l(i)*(x_G(i) - x_end(i-1))/p.d(i);
                y_end(i) = y_end(i-1) + p.l(i)*(y_G(i) - y_end(i-1))/p.d(i);
                eqns_con(i,1) = diff(x_G(i) - x_end(i-1) - p.d(i)*cos(theta(i)),2);
                eqns_con(i,2) = diff(y_G(i) - y_end(i-1) - p.d(i)*sin(theta(i)),2);
            end
        end
        eqns = [reshape(eqns_LMB,[2*p.N 1]); eqns_AMB.'; reshape(eqns_con,[2*p.N 1])];
    end

    % Lagrange equations method
    if strcmp(p.method,'minimal_Lagrange')
        disp('    using minimal coordinates and Lagrange equations...')

        % Initialize symbolic variables
        syms t 
        theta = sym([]);
        theta_p = sym([]);
        dtheta_p = sym([]);
        er = sym([]);
        et = sym([]);
        rG_O = sym([]);
        rG_A = sym([]);
        rA_B = sym([]);
        v_G = sym([]);
        Ek = sym();
        Ep = sym();
        eqns = sym([]);

        % Distance and velocity vectors, kinetic and potential energy
        disp('Creating distance and velocity vectors...')
        for i = 1:p.N
            theta(i) = symfun(sym(sprintf('theta_%d(t)', i)), t);
            theta_p(i) = sym(sprintf('theta_p_%d', i));
            dtheta_p(i) = sym(sprintf('dtheta_p_%d', i));
            er(i,:) = [cos(theta_p(i)) sin(theta_p(i)) 0];
            et(i,:) = cross(k_hat,er(i,:));
            rA_B(i,:) = p.l(i)*er(i,:);
            rG_A(i,:) = p.d(i)*er(i,:);
            rG_O(i,:) = rG_A(i,:);
            if (i > 1)
                v_G(i,:) = v_G(i-1,:) + dtheta_p(i-1)*(p.l(i-1)-p.d(i-1))*et(i-1,:)...
                    + dtheta_p(i)*p.d(i)*et(i,:);
                for j = 1:(i-1)
                    rG_O(i,:) = rG_O(i,:) + rA_B(j,:);
                end
            else
                v_G(1,:) = dtheta_p(1)*p.d(1)*et(1,:);
            end
            Ek = Ek + 1/2*(p.m(i)*sum(v_G(i,:).^2) + p.I_G(i)*dtheta_p(i)^2);
            Ep = Ep - p.m(i)*p.g*rG_O(i,:)*i_hat.';
            if p.springs
                if i == 1
                    Ep = Ep + 1/2*p.k(1)*theta_p(1)^2;
                else
                    Ep = Ep + 1/2*p.k(i)*(theta_p(i) - theta_p(i-1))^2;
                end
            end
        end

        % Lagrange equations
        disp('Writing Lagrange equations...')
        L = Ek - Ep;
        for i = 1:p.N
            eqns(i) = subs(diff(L,theta_p(i)),[theta_p dtheta_p],[theta diff(theta)])...
            - diff(subs((diff(L,dtheta_p(i))),[theta_p dtheta_p],[theta diff(theta)]),t);
        end
    end

    % Integrate the equations in ODE45
    if (strcmp(p.method,'minimal_AMB') || strcmp(p.method,'minimal_Lagrange'))
        vars = theta;
    elseif strcmp(p.method,'maximal')
        vars = [R_x R_y theta x_G y_G];
    end
    [new_eqns,new_vars] = reduceDifferentialOrder(eqns,vars);
    [M,b] = massMatrixForm(new_eqns,new_vars);
    if strcmp(p.method,'maximal')
        M_1 = M(:,2*p.N+1:end);
        c = sym('c%d', [1 2*p.N]);
        A = equationsToMatrix(subs(new_eqns(1:3*p.N),[R_x R_y],c),c);
        M_2 = [A; zeros(numel(new_vars)-numel(A(:,1)),2*p.N)];
        M_new = [M_2 M_1];
        q.all_vars = new_vars(2*p.N+1:end);
        q.M = M_new;
        q.b = subs(b,[R_x R_y],zeros(1,2*p.N));
        q.offset = 2*p.N;
    else
        q.all_vars = new_vars;
        q.M = M;
        q.b = b;
        q.offset = 0;
    end
    disp('Integrating in ODE45...')
    q.num_vars = numel(q.all_vars);
    
    % Functions odeprog.m and odeabort.m created by Tim Franklin
    % Downloaded from MATLAB Central File Exchange
    opts = odeset('OutputFcn',@odeprog,'Events',@odeabort);
    
    if (strcmp(p.method,'minimal_AMB') || strcmp(p.method,'minimal_Lagrange'))
        q.type = 1;
    elseif strcmp(p.method,'maximal')
        q.type = 0;
        x_G_0 = zeros(1,p.N);
        y_G_0 = zeros(1,p.N);
        v_G_x_0 = zeros(1,p.N);
        v_G_y_0 = zeros(1,p.N);
        for i = 1:p.N
            if (i == 1)
                x_G_0(1) = p.d(1)*cos(p.icv(1));
                y_G_0(1) = p.d(1)*sin(p.icv(1));
                v_G_x_0(1) = -p.d(1)*p.icv(p.N+1)*sin(p.icv(1));
                v_G_y_0(1) = p.d(1)*p.icv(p.N+1)*cos(p.icv(1));
            else
                x_G_0(i) = x_G_0(i-1) + (p.l(i-1) - p.d(i-1))*cos(p.icv(i-1))...
                    + p.d(i)*cos(p.icv(i));
                y_G_0(i) = y_G_0(i-1) + (p.l(i-1) - p.d(i-1))*sin(p.icv(i-1))...
                    + p.d(i)*sin(p.icv(i));
                v_G_x_0(i) = v_G_x_0(i-1)...
                    - (p.l(i) - p.d(i))*p.icv(p.N+i-1)*sin(p.icv(i-1))...
                    - p.d(i)*p.icv(p.N+i)*sin(p.icv(i));
                v_G_y_0(i) = v_G_y_0(i-1)...
                    + (p.l(i) - p.d(i))*p.icv(p.N+i-1)*cos(p.icv(i-1))...
                    + p.d(i)*p.icv(p.N+i)*cos(p.icv(i));
            end
        end
        p.icv = [p.icv(1:p.N) x_G_0 y_G_0 p.icv(p.N+1:2*p.N) v_G_x_0 v_G_y_0].';
    end
    [t,z] = ode45(@p48_RHS,tspan,p.icv,opts,q);

    % Store the variables
    ftheta = zeros(numel(t),p.N);
    fdtheta = zeros(numel(t),p.N);
    fx_G = zeros(numel(t),p.N);
    fy_G = zeros(numel(t),p.N);
    fx_end = zeros(numel(t),p.N);
    fy_end = zeros(numel(t),p.N);
    fer = zeros(numel(t),3,p.N);
    fet = zeros(numel(t),3,p.N);
    fv_G = zeros(numel(t),3,p.N);

    % Angles of links
    for i = 1:p.N
        if (strcmp(p.method,'minimal_AMB') || strcmp(p.method,'minimal_Lagrange'))
            ftheta(:,i) = z(:,i);
            fdtheta(:,i) = z(:,i+p.N);
        elseif strcmp(p.method,'maximal')
            ftheta(:,i) = z(:,i);
            fdtheta(:,i) = z(:,i+3*p.N);
        end
    end

    % Positions and velocities of link ends and centers of mass
    for i = 1:p.N
        fx = 0;
        fy = 0;
        for j = 1:i-1
            fx = fx + p.l(j)*cos(ftheta(:,j));
            fy = fy + p.l(j)*sin(ftheta(:,j));
        end
        fx_end(:,i) = fx + p.l(i)*cos(ftheta(:,i));
        fy_end(:,i) = fy + p.l(i)*sin(ftheta(:,i));
        fer(:,:,i) = [cos(ftheta(:,i)) sin(ftheta(:,i)) zeros(numel(t),1)];
        if strcmp(p.method,'maximal')
            fx_G(:,i) = z(:,i+p.N);
            fy_G(:,i) = z(:,i+2*p.N);
            fv_G(:,:,i) = [z(:,i+4*p.N) z(:,i+5*p.N) zeros(numel(t),1)];
        else
            fx_G(:,i) = fx + p.d(i)*cos(ftheta(:,i));
            fy_G(:,i) = fy + p.d(i)*sin(ftheta(:,i));
            for j = 1:numel(t)
                fet(j,:,i) = cross(k_hat,fer(j,:,i));
                if (i == 1)
                    fv_G(j,:,1) = fdtheta(j,1)*p.d(1)*fet(j,:,1);
                else
                    fv_G(j,:,i) = fv_G(j,:,i-1)...
                        + fdtheta(j,i-1)*(p.l(i-1)-p.d(i-1))*fet(j,:,i-1)...
                        + fdtheta(j,i)*p.d(i)*fet(j,:,i);
                end
            end
        end
    end

    % Output time taken to solve
    time = toc;
    if (time < 60)
        disp(['Took ',num2str(toc),' seconds to solve.'])
    else
        disp(['Took ',num2str(toc/60),' minutes to solve.'])
    end
end