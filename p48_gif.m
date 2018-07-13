% Michael R. Buche
% MAE 5730 - Problem 48
% Final Computation Project
% Make GIF of animation

function p48_gif(t,x_end,y_end,p)

    % Make sure too many points are not selected
    if (p.num_points > numel(t))
        disp('Too many points selected for GIF.')
        return
    end
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

    % Interpolate times and positions to evenly spaced and less points
    t_int = linspace(t(1),t(end),p.num_points);
    delay_time = t_int(3) - t_int(2);
    fx_end_int = zeros(numel(t_int),p.N);
    fy_end_int = zeros(numel(t_int),p.N);
    for i = 1:p.N
        fx_end_int(:,i) = interp1(t,x_end(:,i),t_int).';
        fy_end_int(:,i) = interp1(t,y_end(:,i),t_int).';
    end

    % Function gif.m created by Chad Greene
    % Downloaded from MATLAB Central File Exchange
    pause(1)
    gif(p.gif_file_name,'DelayTime',delay_time,'frame',gcf);

    % Update positions of the links and hinges frame by frame
    for j = 1:numel(t_int)
        for i = 1:p.N
            if i == 1
                set(link(i),'xdata',[0 fy_end_int(j,1)],'ydata',[0 fx_end_int(j,1)]);
            else
                set(hinge(i),'xdata',fy_end_int(j,i-1),'ydata',fx_end_int(j,i-1))
                set(link(i),'xdata',[fy_end_int(j,i-1) fy_end_int(j,i)],...
                    'ydata',[fx_end_int(j,i-1) fx_end_int(j,i)])
            end
        end
        drawnow;
        % Update timer if engaged
        if p.timer
            set(timer_text,'string',sprintf('$t = %.2f$',t_int(j)));
        end
        % Make gif
        gif;
    end
end