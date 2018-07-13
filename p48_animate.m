% Michael R. Buche
% MAE 5730 - Problem 48
% Final Computation Project
% Animate the motion

function p48_animate(t,x_end,y_end,p)   

    % Create figure and plot axes
    figure;
    grid on
    grid minor
    axis equal
    axis([-1 1 -1 1]*1.1*sum(p.l))
    hold on
    set(gca,'YDir','reverse')
    xlabel('$y(t)$','Interpreter','LaTeX');
    ylabel('$x(t)$','Interpreter','LaTeX');
    set(gca,'fontsize',p.fontsize);

    % Create links and hinges
    link = zeros(1,p.N);
    hinge = zeros(1,p.N);
    line_width = 2;
    if p.hinges
        hinge_size = 5;
    else
        hinge_size = 1;
    end
    link(1) = plot([0 y_end(1,1)],[0 x_end(1,1)],'k','linewidth',line_width);
    hinge(1) = plot(0,0,'ko','markersize',hinge_size,...
        'MarkerFaceColor','w','linewidth',line_width);
    for i = 2:p.N
        link(i) = plot([y_end(1,i-1) y_end(1,i)],[x_end(1,i-1) x_end(1,i)],...
            'k','linewidth',line_width);
        hinge(i) = plot(y_end(1,i-1),x_end(1,i-1),'ko','markersize',hinge_size,...
            'MarkerFaceColor','w','linewidth',line_width);
    end

    % Create timer textbox if opted to do so
    if p.timer
        timer_text = annotation(gcf,'textbox',[0.825 0.75 0.1 0.1],...
            'string',{'$t = 0.00$'},'fontsize',p.fontsize,'Interpreter','LaTeX');
    end

    % Update positions of the links and hinges in real time
    pause(1)
    xint = zeros(1,p.N);
    yint = zeros(1,p.N);
    time = 0;
    tic

    while time < t(end)/p.animation_timescale
        for i = 1:p.N
            xint(i) = interp1(t,x_end(:,i),time*p.animation_timescale);
            yint(i) = interp1(t,y_end(:,i),time*p.animation_timescale);
            if i == 1
                set(link(i),'xdata',[0 yint(1)],'ydata',[0 xint(1)]);
            else
                set(hinge(i),'xdata',yint(i-1),'ydata',xint(i-1))
                set(link(i),'xdata',[yint(i-1) yint(i)],'ydata',[xint(i-1) xint(i)])
            end
        end
        drawnow;
        time = toc;
        % Update timer
        if p.timer
            set(timer_text,'string',sprintf('$t = %.2f$',time*p.animation_timescale));
        end
    end
end