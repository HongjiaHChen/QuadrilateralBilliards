% This function draws the trajectory of a billiard in a rectangular table.
% User inputs the L_H, L_V, which are the base and height of rectangle respectively.
% As well as side, alpha, position vectors that is a result of running
% square_billiards_cobweb.m

function P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position)
    eps = 1e-10;  % for floating point comparison. Substitute for ==
    % We need to convert the values of P_n (position) onto a rectangular grid now.
    N = length(position);
    % The rectangle lies in the opposing corners of (0,0) and (L_H, L_V).
    P_anim = [position; zeros(size(position))]';   % matrix of x and y co-ords as columns

    for i = 1:N
        if abs(side(i) - 0) < eps
            P_anim(i, :) = [position(i); 0];   % for efficiency purposes we can actually skip this as the 
                                              % matrix remains the same
        elseif abs(side(i)-L_H)< eps
            P_anim(i, :) = [L_H; position(i) - L_H];

        elseif abs(side(i) - (L_H + L_V)) < eps
            P_anim(i, :) = [-position(i) + 2*L_H + L_V; L_V];

        elseif abs(side(i)-(2*L_H+L_V))< eps
            P_anim(i, :) = [0; -position(i) + 2*L_H + 2*L_V];
        end
    end

    plot(P_anim(:,1), P_anim(:,2), 'k', 'linewidth', 2), hold on, pbaspect(double([L_H L_V 1]))  % aspect ratio
    
    % FOR PURPOSES OF BIFURCATION PLOT. DELETE AND UNCOMMENT ABOVE LATER.
    %plot(P_anim(:,1), P_anim(:,2), 'k', 'linewidth', 2), hold on, pbaspect(double([L_H 1.5 1]))  % aspect ratio
    %
    
    
    % Blue square outline
    plot([0 L_H], [0 0], 'b', 'linewidth', 3) % note we just need the vertices
    plot([L_H L_H], [0 L_V], 'b', 'linewidth', 3)
    plot([0 0], [0 L_V], 'b', 'linewidth', 3)
    plot([0 L_H], [L_V L_V], 'b', 'linewidth', 3)

    % Add label for the first square
    text(-0.15, 0, 'A', 'Color','b', 'FontSize', 22)
    text(L_H + 0.05, 0, 'B', 'Color','b', 'FontSize', 22)
    text(L_H + 0.05, L_V, 'C', 'Color','b', 'FontSize', 22)
    text(-0.15, L_V, 'D', 'Color','b', 'FontSize', 22)

    set(gca,'XLim',double([0 L_H]), 'YLim', double([0 L_V]),'xtick', [], 'ytick', []) 
    set(gca, 'XColor','w', 'YColor','w')
    %set(gca,'XLim',double([0 L_H]), 'YLim', double([0 1.5]),'xtick', [], 'ytick', []) % FOR PLOT
    set(gcf,'color','w');  % background color to white instead of usual grey
    
    % Formatting title
    %th = title(sprintf('Rectangle Dimensions: %0.2f x %0.5f. Initial angle: %0.2f. Initial position of %0.2f', ...
    %                                    L_H, L_V, alpha(1), position(1)));
                                        
    th = title(sprintf('Rectangle Dimensions: %0.2f x %0.5f. Initial angle: %0.2f. Initial position of %0.2f', ...
                                        L_H, L_V, alpha(1), position(1)));
    
    titlePos = get( th , 'position');
    titlePos(2) = L_V + 0.04;
    %titlePos(2) = 1.5 + 0.04;      % FOR BIFN DIAGRAM
    set( th , 'position' , titlePos);

    % text(-0.17, 1.04, '(d)', 'FontSize', 22) %for the dissertation
end
