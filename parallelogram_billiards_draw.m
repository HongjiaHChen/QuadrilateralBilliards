% This function draws the trajectory of a billiard in a rectangular table.
% User inputs the L_H, L_V, which are the base and height of rectangle respectively.
% As well as side, alpha, position vectors that is a result of running
% square_billiards_cobweb.m

function P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position)
    eps = 1e-10;  % for floating point comparison. Substitute for ==
    % We need to convert the values of P_n (position) onto a rectangular grid now.
    N = length(position);
    
    P_anim = [position; zeros(size(position))]';   % matrix of x and y co-ords as columns

    for i = 1:N
        if abs(side(i) - 0) < eps
            P_anim(i, :) = [position(i); 0];   % for efficiency purposes we can actually skip this as the 
                                              % matrix remains the same
        elseif abs(side(i)-1) < eps
            P_anim(i, :) = [1+(position(i)-1)*cos(gamma); (position(i)-1)*sin(gamma)];

        elseif abs(side(i)-(1+h/sin(gamma))) < eps
            P_anim(i, :) = [-(position(i)-1-h/sin(gamma))+1+h/tan(gamma); h];

        elseif abs(side(i)-(2+h/sin(gamma)))< eps
            P_anim(i, :) = [(2+2*h/sin(gamma)-position(i))*cos(gamma); (2+2*h/sin(gamma)-position(i))*sin(gamma)];
        end
    end

    plot(P_anim(:,1), P_anim(:,2), 'k', 'linewidth', 2), hold on,  pbaspect(double([1+h/tan(gamma) h 1]))  % aspect ratio
    
    % FOR PURPOSES OF BIFURCATION PLOT. DELETE AND UNCOMMENT ABOVE LATER.
    %plot(P_anim(:,1), P_anim(:,2), 'k', 'linewidth', 2), hold on, pbaspect(double([L_H 1.5 1]))  % aspect ratio
    %
    
    
    % Blue square outline
%     plot([0 1], [0 0], 'b', 'linewidth', 3) % note we just need the vertices
%     plot([1 1+h/mod(tan(gamma),pi)], [0 h], 'b', 'linewidth', 3)
%     plot([h/mod(tan(gamma),pi) 1+h/mod(tan(gamma),pi)], [h h], 'b', 'linewidth', 3)
%     plot([0 h/mod(tan(gamma),pi)], [0 h], 'b', 'linewidth', 3)

    plot([0 1], [0 0], 'b', 'linewidth', 3) % note we just need the vertices
    plot([1 1+h/tan(gamma)], [0 h], 'b', 'linewidth', 3)
    plot([h/tan(gamma) 1+h/tan(gamma)], [h h], 'b', 'linewidth', 3)
    plot([0 h/tan(gamma)], [0 h], 'b', 'linewidth', 3)

    % Add label for the first square
    text(-0.15, 0, 'A', 'Color','b', 'FontSize', 22)
    text(1 + 0.05, 0, 'B', 'Color','b', 'FontSize', 22)
    text(1+h/tan(gamma) + 0.05, h, 'C', 'Color','b', 'FontSize', 22)
    text(h/tan(gamma)-0.15, h, 'D', 'Color','b', 'FontSize', 22)

    axis equal % Hinke's suggestion
    
    set(gca,'XLim',double([0 1+h/tan(gamma)]), 'YLim', double([0 h]),'xtick', [], 'ytick', []) 
    set(gca, 'XColor','w', 'YColor','w')
    %set(gca,'XLim',double([0 L_H]), 'YLim', double([0 1.5]),'xtick', [], 'ytick', []) % FOR PLOT
    set(gcf,'color','w');  % background color to white instead of usual grey
    
    
    % Formatting title
    %th = title(sprintf('Rectangle Dimensions: %0.2f x %0.5f. Initial angle: %0.2f. Initial position of %0.2f', ...
    %                                    L_H, L_V, alpha(1), position(1)));
                                        
    th = title(sprintf('h = %0.2f , gamma = %0.5f. init angle = %0.5f. init pos =  %0.5f', ...
                                        h, gamma, alpha(1), position(1)));
    
    titlePos = get( th , 'position');
    
    
    h = 100; % DELETE THIS. TEMP
    
    
    titlePos(2) = h + 0.04;
    %titlePos(2) = 1.5 + 0.04;      % FOR BIFN DIAGRAM
    set( th , 'position' , titlePos);

    % text(-0.17, 1.04, '(d)', 'FontSize', 22) %for the dissertation
end
