% Michael R. Buche
% MAE 5730 - Problem 48
% Final Computation Project
% Plot the motion

function p48_plot(t,x_end,y_end,theta,p)

    % Plot the angles of each link if any turn over
    links_wrap_around = max(max(abs(theta - wrapToPi(theta)))) > 0;
    if links_wrap_around
        figure
        for i = 1:p.N
            plot(t,rad2deg(theta(:,i)),'linewidth',1);
            hold on
            leg{i} = ['$$\theta_' sprintf('{%d}',i) '(t)$$'];
        end
        hold off
        legendo = legend(leg);
        set(legendo,'Interpreter','latex','FontSize',p.fontsize)
        xlabel('$t$','Interpreter','LaTeX');
        ylabel('$\theta_i(t)$','Interpreter','LaTeX');
        set(gca,'FontSize',p.fontsize);
    end
    
    % Plot the angles of each link wrapped to pi
    figure
    for i = 1:p.N
        plot(t,rad2deg(wrapToPi(theta(:,i))),'linewidth',1);
        hold on
        leg{i} = ['$$\theta_' sprintf('{%d}',i) '(t)$$'];
    end
    hold off
    legendo = legend(leg);
    set(legendo,'Interpreter','latex','FontSize',p.fontsize)
    xlabel('$t$','Interpreter','LaTeX');
    ylabel('$\theta_i(t)$','Interpreter','LaTeX');
    set(gca,'FontSize',p.fontsize);

    % Plot the path of the ends of each link
    figure
    for i = 1:p.N
        plot(y_end(:,i),x_end(:,i),'linewidth',1);
        hold on
        leg{i} = sprintf('%d',i);
    end
    hold off
    set(gca,'YDir','reverse')
    legendo = legend(leg);
    set(legendo,'Interpreter','latex','FontSize',p.fontsize)
    xlabel('$y(t)$','Interpreter','LaTeX');
    ylabel('$x(t)$','Interpreter','LaTeX');
    set(gca,'FontSize',p.fontsize);
    axis equal

    % Plot the path of the end of the last link
    figure
    plot(y_end(:,end),x_end(:,end),'k','linewidth',1);
    legun = sprintf('%d',p.N);
    set(gca,'YDir','reverse')
    legendo = legend(legun);
    set(legendo,'Interpreter','latex','FontSize',p.fontsize)
    xlabel('$y(t)$','Interpreter','LaTeX');
    ylabel('$x(t)$','Interpreter','LaTeX');
    set(gca,'FontSize',p.fontsize);
    axis equal
end