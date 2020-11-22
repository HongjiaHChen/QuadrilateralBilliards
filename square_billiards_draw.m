% This function draws the trajectory of a billiard in a square table.
% User inputs the side, alpha, position that is a result of running
% square_billiards_cobweb.m

function P_anim = square_billiards_draw(side, alpha, position, animflag)
    % We need to convert the values of P_n (position) onto a square grid now.
    N = length(position);
    P_anim = [position; zeros(size(position))]';   % matrix of x and y co-ords as columns

    for i = 1:N
        if side(i) == 0
            P_anim(i, :) = [position(i); 0];   % for efficiency purposes we can actually skip this as the 
                                              % matrix remains the same
        elseif side(i) == 1
            P_anim(i, :) = [1; position(i) - 1];

        elseif side(i) == 2
            % P_anim(i, :) = [position(i) - 2; 1];  %this is a mistake but
            % leads to very pretty pictures
            P_anim(i, :) = [-position(i) + 3; 1];

        elseif side(i) == 3
            % P_anim(i, :) = [0; position(i) - 3]; such as "Using straight lines to draw
            % parabolic curves" from primary school ( for atan(23/37))

            P_anim(i, :) = [0; -position(i) + 4];
        end
    end
    plot(1,1); hold on % TEMP PLOT
    % Blue square outline
    plot([0 1], [0 0], 'b', 'linewidth', 3) % note we just need the vertices
    plot([1 1], [0 1], 'b', 'linewidth', 3)
    plot([0 0], [0 1], 'b', 'linewidth', 3)
    plot([0 1], [1 1], 'b', 'linewidth', 3)

    % Add label for the first square
    text(-0.1, 0, 'A', 'Color','b', 'FontSize', 22)
    text(1.05, 0, 'B', 'Color','b', 'FontSize', 22)
    text(1.05, 1, 'C', 'Color','b', 'FontSize', 22)
    text(-0.1, 1, 'D', 'Color','b', 'FontSize', 22)

    %set(gca, 'XLim',[0 1], 'YLim', [0 1], 'visible','off')
    pbaspect([1 1 1]) % aspect ratio
    set(gca,'XLim',[0 1], 'YLim', [0 1],'xtick', [], 'ytick', [])

    % Formatting title
    %th = title(sprintf('Trajectory for initial angle of %0.2f and initial position of %0.2f', init_angle, init_pos))
    %titlePos = get( th , 'position');
    %titlePos(2) = 1.05;
    %set( th , 'position' , titlePos);

    set(gcf,'color','w');  % background color to white instead of usual grey
    
    % STATIC IMAGE
    if nargin <= 3
        plot(P_anim(:,1), P_anim(:,2), 'k', 'linewidth', 2), hold on
        %text(-0.17, 1.04, '(c)', 'FontSize', 22) %for the dissertation
    else  % animflag has been provided
    
    % DO THE ANIMATIONS
    h = animatedline('Color', 'k', 'LineWidth', 2);
    % Actual animation of trajectory
    for k = 1:N
        addpoints(h, P_anim(k, 1), P_anim(k, 2));
        % drawnow   % draws too fast
        pause(0.4)
    end
    
end

