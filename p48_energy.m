% Michael R. Buche
% MAE 5730 - Problem 48
% Final Computation Project
% Plot total energy

function p48_energy(t,x_G,v_G,theta,dtheta,p)

    % Kinetic and potential energy
    Ek = zeros(1,numel(t));
    Ep = zeros(1,numel(t));
    for i = 1:numel(t)
       Ek(i) = 1/2*(p.m*squeeze(dot(v_G(i,:,:),v_G(i,:,:))) + p.I_G*dtheta(i,:).'.^2);
       Ep(i) = -p.m*p.g*x_G(i,:).';
       if p.springs
           for j = 1:p.N
               if j == 1
                   Ep(i) = Ep(i) + 1/2*p.k(1)*(theta(i,1))^2;
               else
                   Ep(i) = Ep(i) + 1/2*p.k(j)*(theta(i,j) - theta(i,j-1))^2;
               end
           end
       end
    end

    % Total energy
    Et = Ek + Ep;

    % Plot energy versus time
    figure
    plot(t,Ek-Ek(1),'k',t,Ep-Ep(1),'k--',t,Et-Et(1),'k:','linewidth',1)
    legendo = legend('$\Delta E_k$','$\Delta E_p$','$\Delta E_t$');
    set(legendo,'Interpreter','latex','FontSize',p.fontsize)
    xlabel('$t$','Interpreter','LaTeX');
    ylabel('$\Delta E_i$','Interpreter','LaTeX');
    set(gca,'fontsize',p.fontsize);
end