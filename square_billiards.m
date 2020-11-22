clear all; close all; clc; 

% Cobwebbing

% Our Difference Equation for Square Billiards:
% x(n+1) = f(x(n)). f is too complex to write here. Refer to written notes.

% We test different values inital angles and positions to see if they 
% exhibit periodic orbits.

N = 10000;   % number of iterations/maps (bounces)

%% Period 4 orbit with lines and dots. Trivial angle of pi/4 and position of 0.5. 
% Follow the mapping from each position (dot in image)
[side alpha position] = square_billiards_cobweb(pi/4, 0.5, N, 'o-');

% A closed "graph" indicates periodicity.


%% Period 4 orbit with  dots. Trivial angle of pi/4 and position of 0.5. 
for i=[0.01:0.1:4]   % a range of starting positions with the SAME initial ANGLE
    square_billiards_cobweb(pi/4, i, N, 'o'); 
end

% We can do manual cobwebbing here to search for periodicity. Note that the
% map here x(n+1)=f(x(n)) is much more complex than in the circular case
% which is indicates by the strange parallel lines.


%% 3D Scatter Plot P_n, alpha_n, P_{n+1}.

% This shows that indeed a (alpha_n, P_n) pair maps to a unique P_{n+1}

for i=[0.01:0.1:4]   % a range of starting positions with the SAME initial ANGLE
    [side_i alpha_i position_i] = square_billiards_cobweb(pi/4, i, N, 'o', false); % false = no plot 
    
    scatter3(alpha_i(1:(N-1)),position_i(1:(N-1)), position_i(2:N)), hold on
end
title('Cobwebbing Square Billiards')
xlabel('\alpha_n')
ylabel('P_n')
zlabel('P_{n+1}')

rotate3d on
set(gca,'XLim',[0 pi], 'YLim', [0 4], 'ZLim',[0 4])

%% Optional add a plane to the above plot
% P_{n+1} = P_{n} plane
[x y] = meshgrid(0:0.01:pi, 0:0.1:4); % Generate x = alpha_n and y = P_{n} data
z = y;     % z = P_{n+1}

s = surf(x,y,z,'FaceAlpha',0.5) %Plot the surface
s.EdgeColor = 'none';


%% Non-periodic orbit with dots and lines. Initial angle of 0.15 and position of 0.5.
square_billiards_cobweb(0.15, 0.5, N, 'o-');
% Notice that we don't have a closed shape. If we still aren't sure,
% increase N. If the picture changes, we have a non-periodic orbit. If it
% reamins the same, then it may be periodic with a very large period.

%% Increasing the number of iterations from the example above
square_billiards_cobweb(0.15, 0.5, 10 * N, 'o-');

% We see that the small 'holes' in the bottom left corner are now filled.
% It is actually difficult to examine if we have periodicity here.

%% Non-periodic orbit with dots. Initial angle of 0.15 and position of 0.5.
[side alpha position] = square_billiards_cobweb(0.15, 0.5, N, 'o');

% Notice the apparent duplication in our mapping. For each P_n we seem to
% have two P_n+1 possible. This is beacause we haven't taken into account
% the angle on this iteration.

%% Plot for the angle
plot(alpha(1:(N-1)), alpha(2:N), 'o'), hold on 
plot([0:0.1:pi], [0:0.1:pi], 'r-')
xlabel('\alpha_n')
ylabel('\alpha_{n+1}')

xticks([0 pi/2 pi])
xticklabels({'0','\pi/2', '\pi'})
yticks([0 pi/2 pi])
yticklabels({'0','\pi/2', '\pi'})

set(gca,'XLim',[0 pi], 'YLim', [0 pi])

%% For Hinke

% Having an angle of 0.15 off all 4 of our sides
for i=[0.01:0.1:4]   % a range of starting positions with the SAME initial ANGLE AND ONE ITERATION
    square_billiards_cobweb(0.15, i, 2, 'o'); 
end

% Creating a grid
% Vertical lines
for i=1:3
    plot(repelem(i, 101), [0:0.04:4], 'Color', [0.7 0.7 0.7], 'LineStyle', '--'), hold on
end
% Horizontal lines
for i=1:3
    plot([0:0.04:4], repelem(i, 101), 'Color', [0.7 0.7 0.7], 'LineStyle', '--'), hold on
end



%% 
[side alpha position] = square_billiards_cobweb(0.15, 0.5, N, 'o');
square_billiards_cobweb(0.15, 0.5, 6, 'o-');

% Our mapping is indeed correct and cobwebbing is indeed working. We just
% aren't sure which line to map to without alpha_n

%% 
[side alpha position] = square_billiards_cobweb(0.15, 0.5, N, 'o');
square_billiards_cobweb(0.15, 0.5, 6, 'o-');
% Creating a grid
% Vertical lines
for i=1:3
    plot(repelem(i, 101), [0:0.04:4], 'Color', [0.7 0.7 0.7], 'LineStyle', '--'), hold on
end
% Horizontal lines
for i=1:3
    plot([0:0.04:4], repelem(i, 101), 'Color', [0.7 0.7 0.7], 'LineStyle', '--'), hold on
end

% Notice that the 'sharp' edges occur at exactly 0,1,2,3. See
% Square_Chaotic_Cobweb.PNG for more information.

%% 3D Scatter Plot
scatter3(alpha(1:(N-1)),position(1:(N-1)), position(2:N))
title('Cobwebbing Square Billiards')
xlabel('\alpha_n')
ylabel('P_n')
zlabel('P_{n+1}')

set(gca,'XLim',[0 pi], 'YLim', [0 4], 'ZLim',[0 4])

xticks([0 pi/2 pi])
xticklabels({'0','\pi/2', '\pi'})

rotate3d on

%% Another non-periodic orbit for the angle of 1.5
[side_chaos alpha_chaos position_chaos] = square_billiards_cobweb(1.5, 0.5, N, 'o');
hold off
%% 
scatter3(alpha_chaos(1:(N-1)),position_chaos(1:(N-1)), position_chaos(2:N))

xlabel('\alpha_n')
ylabel('P_n')
zlabel('P_{n+1}')

set(gca,'XLim',[0 pi], 'YLim', [0 4], 'ZLim',[0 4])

xticks([0 pi/2 pi])
xticklabels({'0','\pi/2', '\pi'})

rotate3d on



%% Period 8 orbit with dots and lines. Initial angle of atan(1/3) and position of 0.5.
square_billiards_cobweb(atan(1/3), 0.5, N, 'o-');
% The "graph" is already very complicated.
% We have a ratio of 1 top/bottom hit to  3 left/right hits. Check notes on
% square billiards and unfolding for the reason.


%% Period 8 orbit with  dots. Angle of atan(1/3) and position of 0.5. 
for i=[0.01:0.1:1]   % a range of starting positions with the SAME initial ANGLE on side 0
    square_billiards_cobweb(atan(1/3), i, N, 'o'); 
end

% Creating a grid
% Vertical lines
for i=1:3
    plot(repelem(i, 101), [0:0.04:4], 'Color', [0.7 0.7 0.7], 'LineStyle', '--'), hold on
end
% Horizontal lines
for i=1:3
    plot([0:0.04:4], repelem(i, 101), 'Color', [0.7 0.7 0.7], 'LineStyle', '--'), hold on
end


%%

% Further example: 
% square_billiards_cobweb(atan(3/2), 0.5, N, 'o-');
% Ratio of 3 top/bottom hit to  2 left/right hits.


%% Some pretty pictures
square_billiards_cobweb(pi/4-0.01, 0.5, N, 'o-');
%%
square_billiards_cobweb(0.001+0.1*2.6, 0.5, 1000, 'o-');
%% 
square_billiards_cobweb(pi/2+0.01, 0.5, N, 'o-');
%%
square_billiards_cobweb(pi/2+0.1, 0.5, N, 'o-');
%%
square_billiards_cobweb(atan(4/7)+0.1, 0.5, 1000, 'o-');
%%
[side alpha position] = square_billiards_cobweb(atan(3/5), 0.5, N, 'o-');







%% 
%% 





%% 



%% 3 DIMENSIONAL P_n, alpha_n and P_{n+1}
%% 
epsilon = 0.001;
[alpha_0 P_0] = meshgrid(0 + epsilon:0.1:pi-epsilon, 0 + epsilon:0.1:4-epsilon); 
% Generate x = alpha_n and y = P_{n} data
% Not good to include starting angles of exactly 0 and pi as these result
% in discontinuties in our plot.

alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

for row=1:size(alpha_0,1) % over the rows
    
    for col=1:size(alpha_0,2) % over the columns
        [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 1, ...
                                                                                    'o', false); % no plot
        alpha_1(row, col) = alpha_temp(2);
        P_1(row, col) = pos_temp(2);
       
    end
end

alpha_1;
P_1;

surf(alpha_0, P_0, P_1) %Plot the surface
xlabel('\alpha_n')
ylabel('P_n')
zlabel('P_{n+1}')

set(gca,'XLim',[0 pi], 'YLim', [0 4], 'ZLim',[0 4])
xticks([0 pi/2 pi])
xticklabels({'0','\pi/2', '\pi'})

rotate on  % For this one, we have an error and click on the rotate button manually


%% Trying to have 4 seperate plots by plotting them seperately and using hold on
epsilon = 0.001;

for side_id = 0:3

    [alpha_0 P_0] = meshgrid(0 + epsilon:0.01:pi-epsilon, side_id + epsilon:0.01:side_id+1-epsilon); 
    % Generate x = alpha_n and y = P_{n} data
    % Not good to include starting angles of exactly 0 and pi as these result
    % in discontinuties in our plot.

    alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

    for row=1:size(alpha_0,1) % over the rows

        for col=1:size(alpha_0,2) % over the columns
            [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 1, ...
                                                                                        'o', false); % no plot
            alpha_1(row, col) = alpha_temp(2);
            P_1(row, col) = pos_temp(2);
        end
    end

    surf(alpha_0, P_0, P_1, 'Linestyle', 'None') %Plot the surface
    xlabel('\alpha_n')
    ylabel('P_n')
    zlabel('P_{n+1}')
    hold on
end

set(gca,'XLim',[0 pi], 'YLim', [0 4], 'ZLim',[0 4])
xticks([0 pi/2 pi])
xticklabels({'0','\pi/2', '\pi'})


% Figuring out the exact locations of the discontinuties that occur when
% P_n is one sides [1,2) and [2,3)

% The disconituites occur when we map from the sides [1,2) or [2,3) to
% corner 4 (or 0 since 4 mod 4 = 0)

% Side [1,2)
discont_P_12 = [1+epsilon:0.01:2-epsilon];
discont_alpha_12 = atan(1./(1-discont_P_12)) + pi;
scatter3(discont_alpha_12, discont_P_12, repelem(4, size(discont_P_12,2)), 'green'), hold on

% Side [2,3)
discont_P_23 = [2+epsilon:0.01:3-epsilon];
discont_alpha_23 = acot(3-discont_P_23);
scatter3(discont_alpha_23, discont_P_23, repelem(4, size(discont_P_23,2)), 'red')

rotate on  % For this one, we have an error and click on the rotate button manually


% Notice that all the plots for each side are parallel and are only shifted
% vertically. To understand what happends, we only need to study a single
% side, P_n in (0,1).



%% Trying to remove the last 2 wall discontinitues by finding the equation for the curve
% When our map eventually hits 4 and goes back down to 0 due to the mod 4.

% For sides 1 and 2 that hit the vertex 0.

% For side 1 -> 0, we require that P_n+1 = P_n + 1  (mod 4).
% Futhermore, the extra "stuff" is 0.
% i.e: 1+x_n*tan(alpha) = 0 => alpha = atan(1/(1-P_n)), 
% where x_n = P_n - floor(P_n), and floor(P_n) = 1.
% We add an additional pi to keep alpha in the correct range as 0<alpha<pi.

% For side 2 -> 0, we require that P_n+1 = P_n + 2  (mod 4).
% Futhermore, the extra "stuff" is 0.
% i.e: 1-x_n-cot(alpha)=0 => alpha = acot(3-P_n), 
% where x_n = P_n - floor(P_n), and floor(P_n) = 2.


side_id = 2;

%discont_P_12 = [1+epsilon:0.01:2-epsilon];
%discont_alpha_12 = atan(1./(1-discont_P_12)) + pi;

[alpha_0 P_0] = meshgrid(0 + epsilon:0.05:pi-epsilon, side_id + epsilon:0.01:side_id+1-epsilon); 
%[alpha_0 P_0] = meshgrid(discont_alpha_12, discont_P_12)

% Generate x = alpha_n and y = P_{n} data
% Not good to include starting angles of exactly 0 and pi as these result
% in discontinuties in our plot.

alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));


% HALF OF THE PLOT, TO AVOID DISCONTINUITY "WALL"

for row=1:size(alpha_0,1) % over the rows

    for col=1:size(alpha_0,2) % over the columns
        
        %discont_alpha_12 = atan(1./(1-P_0(row, col))) + pi;
        discont_alpha_23 = acot(3-P_0(row,col));
        
        if alpha_0(row,col) < discont_alpha_23             % just the first half of the plot

            [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 1, ...
                                                                                        'o', false); % no plot
            alpha_1(row, col) = alpha_temp(2);
            P_1(row, col) = pos_temp(2);
        else
            alpha_1(row, col) = NaN;
            P_1(row, col) = NaN;
            
        end
    end
end

%surf(alpha_0, P_0, P_1, 'Linestyle', 'None') %Plot the surface
surf(alpha_0, P_0, P_1) %Plot the surface
hold on

% OTHER HALF OF THE PLOT

for row=1:size(alpha_0,1) % over the rows
    for col=1:size(alpha_0,2) % over the columns
        
        % discont_alpha_12 = atan(1./(1-P_0(row, col))) + pi;
        discont_alpha_23 = acot(3-P_0(row, col));
        
        if alpha_0(row,col) >= discont_alpha_23             % second half of the plot

            [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 1, ...
                                                                                        'o', false); % no plot
            alpha_1(row, col) = alpha_temp(2);
            P_1(row, col) = pos_temp(2);
        else
            alpha_1(row, col) = NaN;
            P_1(row, col) = NaN;
            
        end
    end
end
xlabel('\alpha_n')
ylabel('P_n')
zlabel('P_{n+1}')
surf(alpha_0, P_0, P_1) %Plot the surface

%% Complete Surface with no wall discontinuities

epsilon = 0.001;

for side_id = 0:3

    [alpha_0 P_0] = meshgrid(0 + epsilon:0.05:pi-epsilon, side_id + epsilon:0.01:side_id+1-epsilon); 
    % Generate x = alpha_n and y = P_{n} data
    % Not good to include starting angles of exactly 0 and pi as these result
    % in discontinuties in our plot.

    alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

    for row=1:size(alpha_0,1) % over the rows
        for col=1:size(alpha_0,2) % over the columns
            [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 2, ...
                                                                                        'o', false); % no plot
                                                                                    
            ith_iterate = 2;   %depends on the N we have chosen , x_1, x_2, x_3,...      
            % change this for any value from 2 to N+1, for period
            % ith_iterate - 1
            
            alpha_1(row, col) = alpha_temp(ith_iterate);
            P_1(row, col) = pos_temp(ith_iterate);
            
            % Fixing discontinuities on side 1 and 2
            % Side 1
            if side_id == 1
                discont_alpha_12 = atan(1./(1-P_0(row, col))) + pi;
                if alpha_0(row,col) >= discont_alpha_12             % second half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN; %ignore these for now, so MATLAB doesn't connect them
                end
                
            elseif side_id == 2
                discont_alpha_23 = acot(3-P_0(row, col));
                if alpha_0(row,col) >= discont_alpha_23             % second half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN; %ignore these for now, so MATLAB doesn't connect them
                    %disp('A')
                end
            end
        end
    end
    surf(alpha_0, P_0, P_1), hold on %Plot the surface
end

% Adding additional surfaces for sides 1 & 2 for which we've omitted to
% avoid the discontiunuity
for side_id = 1:2
    [alpha_0 P_0] = meshgrid(0 + epsilon:0.05:pi-epsilon, side_id + epsilon:0.01:side_id+1-epsilon); 
    alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

    for row=1:size(alpha_0,1) % over the rows
        for col=1:size(alpha_0,2) % over the columns

            [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 1, ...
                                                                                            'o', false); % no plot
            alpha_1(row, col) = alpha_temp(2);
            P_1(row, col) = pos_temp(2);
            
            if side_id == 1
                discont_alpha_12 = atan(1./(1-P_0(row, col))) + pi;
                if alpha_0(row,col) < discont_alpha_12             % omit first half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN;  % we already drew these
                end
                
            elseif side_id == 2
                discont_alpha_23 = acot(3-P_0(row, col));
                if alpha_0(row,col) < discont_alpha_23             % omit first half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN; % we already drew these
                end
            end
        end
    end
    surf(alpha_0, P_0, P_1) %Plot the surface
end



set(gca,'XLim',[0 pi], 'YLim', [0 4], 'ZLim',[0 4])
xticks([0 pi/2 pi])
xticklabels({'0','\pi/2', '\pi'})
%title(sprintf('Period: %d', ith_iterate-1))




%% Just the single side, P_n in (0,1)

epsilon = 0.001;
[alpha_0 P_0] = meshgrid(0 + epsilon:0.01:pi-epsilon, 0 + epsilon:0.1:1-epsilon); 
% Generate x = alpha_n and y = P_{n} data
% Not good to include starting angles of exactly 0 and pi as these result
% in discontinuties in our plot.

alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

for row=1:size(alpha_0,1) % over the rows
    
    for col=1:size(alpha_0,2) % over the columns
        [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 1, ...
                                                                                    'o', false); % no plot
        alpha_1(row, col) = alpha_temp(2);
        P_1(row, col) = pos_temp(2);
       
    end
end

alpha_1;
P_1;

surf(alpha_0, P_0, P_1) %Plot the surface
xlabel('\alpha_n')
ylabel('P_n')
zlabel('P_{n+1}')

set(gca,'XLim',[0 pi], 'YLim', [0 1], 'ZLim',[0 4])
xticks([0 pi/2 pi])
xticklabels({'0','\pi/2', '\pi'})
title('Just side 0')
rotate on  % For this one, we have an error and click on the rotate button manually

%% TEMP FOR RElATIONSHIP

% Cobweb diagram atan(0.15)
N_temp = 1500;
[side alpha position] = square_billiards_cobweb(0.3, 0.5, N_temp, 'o', false);
scatter3(alpha(1:(N_temp-1)),position(1:(N_temp-1)), position(2:N_temp), 'black')
%line(alpha(1:(N_temp-1)), position(1:(N_temp-1)), position(2:N_temp))
hold on

xlabel('\alpha_n')
ylabel('P_n')
zlabel('P_{n+1}')

set(gca,'XLim',[0 pi], 'YLim', [0 4], 'ZLim',[0 4])
xticks([0 pi/2 pi])
xticklabels({'0','\pi/2', '\pi'})

%plot()
% Cobweb surface

epsilon = 0.001;

for side_id = 0:3

    [alpha_0 P_0] = meshgrid(0 + epsilon:0.05:pi-epsilon, side_id + epsilon:0.01:side_id+1-epsilon); 
    % Generate x = alpha_n and y = P_{n} data
    % Not good to include starting angles of exactly 0 and pi as these result
    % in discontinuties in our plot.

    alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

    for row=1:size(alpha_0,1) % over the rows
        for col=1:size(alpha_0,2) % over the columns
            [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 2, ...
                                                                                        'o', false); % no plot
                                                                                    
            ith_iterate = 2;   %depends on the N we have chosen , x_1, x_2, x_3,...      
            % change this for any value from 2 to N+1, for period
            % ith_iterate - 1
            
            alpha_1(row, col) = alpha_temp(ith_iterate);
            P_1(row, col) = pos_temp(ith_iterate);
            
            % Fixing discontinuities on side 1 and 2
            % Side 1
            if side_id == 1
                discont_alpha_12 = atan(1./(1-P_0(row, col))) + pi;
                if alpha_0(row,col) >= discont_alpha_12             % second half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN; %ignore these for now, so MATLAB doesn't connect them
                end
                
            elseif side_id == 2
                discont_alpha_23 = acot(3-P_0(row, col));
                if alpha_0(row,col) >= discont_alpha_23             % second half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN; %ignore these for now, so MATLAB doesn't connect them
                    %disp('A')
                end
            end
        end
    end
    %surf(alpha_0, P_0, P_1, 'Linestyle', 'None'), hold on %Plot the surface
    surf(alpha_0, P_0, P_1, alpha_1, 'Linestyle', 'None'), hold on % No mesh on surface
    % color depends on alpha_1
    
    cmap = [ cool(64); winter(64) ];
    colormap(gca, cmap);
    clims = [0 pi];  % minimum and maximum alpha values
    caxis(clims);
    
    c = colorbar('Limits', [0 pi], 'Ticks', [0 pi], 'TickLabels', {'0', '\pi'}, 'FontSize', 10)
    c.Label.String = '\alpha_{n+1}'; c.Label.FontSize = 15;
    c.Label.Rotation = 0; % to rotate the label
    pos = get(c,'Position');
    c.Label.Position = [pos(1)+1.5 pos(2)+pi/2]; % to change its position
end

% Adding additional surfaces for sides 1 & 2 for which we've omitted to
% avoid the discontiunuity
for side_id = 1:2
    [alpha_0 P_0] = meshgrid(0 + epsilon:0.05:pi-epsilon, side_id + epsilon:0.01:side_id+1-epsilon); 
    alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

    for row=1:size(alpha_0,1) % over the rows
        for col=1:size(alpha_0,2) % over the columns

            [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 1, ...
                                                                                            'o', false); % no plot
            alpha_1(row, col) = alpha_temp(2);
            P_1(row, col) = pos_temp(2);
            
            if side_id == 1
                discont_alpha_12 = atan(1./(1-P_0(row, col))) + pi;
                if alpha_0(row,col) < discont_alpha_12             % omit first half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN;  % we already drew these
                end
                
            elseif side_id == 2
                discont_alpha_23 = acot(3-P_0(row, col));
                if alpha_0(row,col) < discont_alpha_23             % omit first half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN; % we already drew these
                end
            end
        end
    end
    %surf(alpha_0, P_0, P_1, 'Linestyle', 'None') %Plot the surface
    surf(alpha_0, P_0, P_1, alpha_1, 'Linestyle', 'None'), hold on % No mesh on surface
    % color depends on alpha_1
end








%% Checking if our cobweb diagram atan(0.15) lies on our cobweb surface

% Cobweb diagram atan(0.15)
[side alpha position] = square_billiards_cobweb(0.15, 0.5, N, 'o', false);
scatter3(alpha(1:(N-1)),position(1:(N-1)), position(2:N), 'red')
hold on


% Another cobweb in 3D with slope atan(1.5)
[side_chaos alpha_chaos position_chaos] = square_billiards_cobweb(1.5, 0.5, N, 'o', false);
scatter3(alpha_chaos(1:(N-1)),position_chaos(1:(N-1)), position_chaos(2:N), 'yellow')
hold on

% Cobweb surface

epsilon = 0.001;

for side_id = 0:3

    [alpha_0 P_0] = meshgrid(0 + epsilon:0.05:pi-epsilon, side_id + epsilon:0.01:side_id+1-epsilon); 
    % Generate x = alpha_n and y = P_{n} data
    % Not good to include starting angles of exactly 0 and pi as these result
    % in discontinuties in our plot.

    alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

    for row=1:size(alpha_0,1) % over the rows
        for col=1:size(alpha_0,2) % over the columns
            [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 2, ...
                                                                                        'o', false); % no plot
                                                                                    
            ith_iterate = 2;   %depends on the N we have chosen , x_1, x_2, x_3,...      
            % change this for any value from 2 to N+1, for period
            % ith_iterate - 1
            
            alpha_1(row, col) = alpha_temp(ith_iterate);
            P_1(row, col) = pos_temp(ith_iterate);
            
            % Fixing discontinuities on side 1 and 2
            % Side 1
            if side_id == 1
                discont_alpha_12 = atan(1./(1-P_0(row, col))) + pi;
                if alpha_0(row,col) >= discont_alpha_12             % second half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN; %ignore these for now, so MATLAB doesn't connect them
                end
                
            elseif side_id == 2
                discont_alpha_23 = acot(3-P_0(row, col));
                if alpha_0(row,col) >= discont_alpha_23             % second half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN; %ignore these for now, so MATLAB doesn't connect them
                    %disp('A')
                end
            end
        end
    end
    %surf(alpha_0, P_0, P_1), hold on %Plot the surface
    surf(alpha_0, P_0, P_1, 'Linestyle', 'None')
    cmap = [ cool(64); winter(64) ];
    colormap(gca, cmap);
    clims = [0 pi];  % minimum and maximum alpha values
    caxis(clims);
end

% Adding additional surfaces for sides 1 & 2 for which we've omitted to
% avoid the discontiunuity
for side_id = 1:2
    [alpha_0 P_0] = meshgrid(0 + epsilon:0.05:pi-epsilon, side_id + epsilon:0.01:side_id+1-epsilon); 
    alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

    for row=1:size(alpha_0,1) % over the rows
        for col=1:size(alpha_0,2) % over the columns

            [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 1, ...
                                                                                            'o', false); % no plot
            alpha_1(row, col) = alpha_temp(2);
            P_1(row, col) = pos_temp(2);
            
            if side_id == 1
                discont_alpha_12 = atan(1./(1-P_0(row, col))) + pi;
                if alpha_0(row,col) < discont_alpha_12             % omit first half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN;  % we already drew these
                end
                
            elseif side_id == 2
                discont_alpha_23 = acot(3-P_0(row, col));
                if alpha_0(row,col) < discont_alpha_23             % omit first half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN; % we already drew these
                end
            end
        end
    end
    %surf(alpha_0, P_0, P_1) %Plot the surface
end

xlabel('\alpha_n')
ylabel('P_n')
zlabel('P_{n+1}')

set(gca,'XLim',[0 pi], 'YLim', [0 4], 'ZLim',[0 4])
xticks([0 pi/2 pi])
xticklabels({'0','\pi/2', '\pi'})
%title(sprintf('Period: %d', ith_iterate-1))

%% alpha_{n+1} Surface

epsilon = 0.001;

for side_id = 0:3

    [alpha_0 P_0] = meshgrid(0 + epsilon:0.05:pi-epsilon, side_id + epsilon:0.01:side_id+1-epsilon); 
    % Generate x = alpha_n and y = P_{n} data
    % Not good to include starting angles of exactly 0 and pi as these result
    % in discontinuties in our plot.

    alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

    for row=1:size(alpha_0,1) % over the rows
        for col=1:size(alpha_0,2) % over the columns
            [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 2, ...
                                                                                        'o', false); % no plot
                                                                                    
            ith_iterate = 2;   %depends on the N we have chosen , x_1, x_2, x_3,...      
            % change this for any value from 2 to N+1, for period
            % ith_iterate - 1
            
            alpha_1(row, col) = alpha_temp(ith_iterate);
            P_1(row, col) = pos_temp(ith_iterate);
            
            % Fixing discontinuities on side 1 and 2
            % Side 1
            if side_id == 1
                discont_alpha_12 = atan(1./(1-P_0(row, col))) + pi;
                if alpha_0(row,col) >= discont_alpha_12             % second half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN; %ignore these for now, so MATLAB doesn't connect them
                end
                
            elseif side_id == 2
                discont_alpha_23 = acot(3-P_0(row, col));
                if alpha_0(row,col) >= discont_alpha_23             % second half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN; %ignore these for now, so MATLAB doesn't connect them
                    %disp('A')
                end
            end
        end
    end
    %surf(alpha_0, P_0, alpha_1), hold on %Plot the surface
    surf(alpha_0, P_0, alpha_1, 'Linestyle', 'None')
end

% Adding additional surfaces for sides 1 & 2 for which we've omitted to
% avoid the discontiunuity
for side_id = 1:2
    [alpha_0 P_0] = meshgrid(0 + epsilon:0.05:pi-epsilon, side_id + epsilon:0.01:side_id+1-epsilon); 
    alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

    for row=1:size(alpha_0,1) % over the rows
        for col=1:size(alpha_0,2) % over the columns

            [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 1, ...
                                                                                            'o', false); % no plot
            alpha_1(row, col) = alpha_temp(2);
            P_1(row, col) = pos_temp(2);
            
            if side_id == 1
                discont_alpha_12 = atan(1./(1-P_0(row, col))) + pi;
                if alpha_0(row,col) < discont_alpha_12             % omit first half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN;  % we already drew these
                end
                
            elseif side_id == 2
                discont_alpha_23 = acot(3-P_0(row, col));
                if alpha_0(row,col) < discont_alpha_23             % omit first half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN; % we already drew these
                end
            end
        end
    end
    %surf(alpha_0, P_0, alpha_1) %Plot the surface
    surf(alpha_0, P_0, alpha_1, 'Linestyle', 'None')
end

xlabel('\alpha_n')
ylabel('P_n')
zlabel('\alpha_{n+1}')

set(gca,'XLim',[0 pi], 'YLim', [0 4], 'ZLim',[0 pi])
xticks([0 pi/2 pi])
xticklabels({'0','\pi/2', '\pi'})
%title(sprintf('Period: %d', ith_iterate-1))

% Note the sharp "triangular" edges, this is again due to us animing at the
% corner. Since the formula for mapping changes dramatically depending on
% the side.


% SOMETHING WRONG. PROB DUE TO THE METHOD USED TO AVOID DISCOUNTINUITY IS
% SPECFIICALLY FOR P_n not Alpha_n


%% P_{n+1} Surface for Dissertation

% We also make the color of the surface correspond to the alpha_{n+1}
% value. Our surface then contains information about all 4 variables.

epsilon = 0.001;

for side_id = 0:3

    [alpha_0 P_0] = meshgrid(0 + epsilon:0.005:pi-epsilon, side_id + epsilon:0.001:side_id+1-epsilon); 
    % Generate x = alpha_n and y = P_{n} data
    % Not good to include starting angles of exactly 0 and pi as these result
    % in discontinuties in our plot.

    alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

    for row=1:size(alpha_0,1) % over the rows
        for col=1:size(alpha_0,2) % over the columns
            [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 2, ...
                                                                                        'o', false); % no plot
                                                                                    
            ith_iterate = 2;   %depends on the N we have chosen , x_1, x_2, x_3,...      
            % change this for any value from 2 to N+1, for period
            % ith_iterate - 1
            
            alpha_1(row, col) = alpha_temp(ith_iterate);
            P_1(row, col) = pos_temp(ith_iterate);
            
            % Fixing discontinuities on side 1 and 2
            % Side 1
            if side_id == 1
                discont_alpha_12 = atan(1./(1-P_0(row, col))) + pi;
                if alpha_0(row,col) >= discont_alpha_12             % second half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN; %ignore these for now, so MATLAB doesn't connect them
                end
                
            elseif side_id == 2
                discont_alpha_23 = acot(3-P_0(row, col));
                if alpha_0(row,col) >= discont_alpha_23             % second half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN; %ignore these for now, so MATLAB doesn't connect them
                    %disp('A')
                end
            end
        end
    end
    surf(alpha_0, P_0, P_1, alpha_1, 'Linestyle', 'None'), hold on % No mesh on surface
    % color depends on alpha_1
    cmap = [ cool(64); winter(64) ];
    colormap(gca, cmap);
    clims = [0 pi];  % minimum and maximum alpha values
    caxis(clims);
    
    c = colorbar('Limits', [0 pi], 'Ticks', [0 pi], 'TickLabels', {'0', '\pi'}, 'FontSize', 10)
    c.Label.String = '\alpha_{n+1}'; c.Label.FontSize = 15;
    c.Label.Rotation = 0; % to rotate the label
    pos = get(c,'Position');
    c.Label.Position = [pos(1)+1.5 pos(2)+pi/2]; % to change its position
end

% Adding additional surfaces for sides 1 & 2 for which we've omitted to
% avoid the discontiunuity
for side_id = 1:2
    [alpha_0 P_0] = meshgrid(0 + epsilon:0.005:pi-epsilon, side_id + epsilon:0.001:side_id+1-epsilon); 
    alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

    for row=1:size(alpha_0,1) % over the rows
        for col=1:size(alpha_0,2) % over the columns

            [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 1, ...
                                                                                            'o', false); % no plot
            alpha_1(row, col) = alpha_temp(2);
            P_1(row, col) = pos_temp(2);
            
            if side_id == 1
                discont_alpha_12 = atan(1./(1-P_0(row, col))) + pi;
                if alpha_0(row,col) < discont_alpha_12             % omit first half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN;  % we already drew these
                end
                
            elseif side_id == 2
                discont_alpha_23 = acot(3-P_0(row, col));
                if alpha_0(row,col) < discont_alpha_23             % omit first half of the plot
                    alpha_1(row, col) = NaN;
                    P_1(row, col) = NaN; % we already drew these
                end
            end
        end
    end
    surf(alpha_0, P_0, P_1, 'Linestyle', 'None') %Plot the surface
    
end

xlabel('\alpha_n', 'FontSize', 15)
ylabel('P_n', 'FontSize', 15)
zlabel('P_{n+1}', 'FontSize', 15)

set(gca,'XLim',[0 pi], 'YLim', [0 4], 'ZLim',[0 4])
xticks([0 pi/2 pi])
xticklabels({'0','\pi/2', '\pi'})

view(-72,24)

%% Investigating higher order iterate of P. Eg: P_n+2, P_n+3, etc.

% Our surface

epsilon = 0.001;

for side_id = 0:3

    [alpha_0 P_0] = meshgrid(0 + epsilon:0.05:pi-epsilon, side_id + epsilon:0.01:side_id+1-epsilon); 
    % Generate x = alpha_n and y = P_{n} data
    % Not good to include starting angles of exactly 0 and pi as these result
    % in discontinuties in our plot.

    alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

    for row=1:size(alpha_0,1) % over the rows

        for col=1:size(alpha_0,2) % over the columns
            [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 5, ...
                                                                                        'o', false); % no plot
                                                                                    
            ith_iterate = 4;   %depends on the N we have chosen , x_1, x_2, x_3,...      
            % change this for any value from 2 to N+1, 
            % for a period  of ith_iterate - 1
            
            alpha_1(row, col) = alpha_temp(ith_iterate);
            P_1(row, col) = pos_temp(ith_iterate);
        end
    end
    
    % Converting into vector as scatter3 wants vectors not matrices
    alpha_0_vec = alpha_0(:); P_0_vec = P_0(:); P_1_vec = P_1(:);
    
    S = repmat([20],numel(alpha_0_vec),1);
    C = repmat([3],numel(alpha_0_vec),1);
    s = S(:);
    c = C(:);

    scatter3(alpha_0_vec, P_0_vec, P_1_vec, c, s)
    %surf(alpha_0, P_0, P_1, 'Linestyle', 'None') %Plot the surface
    %surf(alpha_0, P_0, P_1) %Plot the surface
    xlabel('\alpha_n')
    ylabel('P_n')
    zlabel(sprintf('P_{n+%d}', ith_iterate-1))
    hold on
end

set(gca,'XLim',[0 pi], 'YLim', [0 1], 'ZLim',[0 4])
xticks([0 pi/2 pi])
xticklabels({'0','\pi/2', '\pi'})
title(sprintf('Period: %d', ith_iterate-1))

%[alpha_0 P_0] = meshgrid(0 + epsilon:0.05:pi-epsilon, 0 + epsilon:0.01:4+1-epsilon); 
%P_1 = P_0;
%plane = surf(alpha_0, P_0, P_1, 'FaceAlpha', 1, 'FaceColor', 'r')


% Our Plane 
[x y] = meshgrid(0:0.01:pi, 0:0.1:4); % Generate x = alpha_n and y = P_{n} data
z = y;     % z = P_{n+1}

s = surf(x,y,z,'FaceAlpha',0.5, 'FaceColor', 'r') %Plot the surface
s.EdgeColor = 'none';

%set(plane,'facealpha',0.2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Investigating period 3 orbits

% A candidate (this initial angle and position are chosen for us to narrowly miss corner 2)
Niter = 3;

[side_3 alpha_3 position_3] = square_billiards_cobweb(pi/4+0.001, 0.001, Niter, 'o-', false);

alpha_3; 
position_3
scatter3(alpha_3(1), position_3(1), position_3(1+Niter), 'blue')

 % Another candidate (this initial angle and position are chosen for us to narrowly miss corner 2)
 [side_3 alpha_3 position_3] = square_billiards_cobweb(acot(1/2)+0.00001, 0.5, Niter, 'o-', false);
 position_3
 scatter3(alpha_3(1), position_3(1), position_3(1+Niter), 'blue')
 
% Another candidate on other surface
% (this initial angle and position are chosen for us to narrowly miss corner 3)
[side_3 alpha_3 position_3] = square_billiards_cobweb(pi - acot(1/2)+0.00001, 0.5, Niter, 'o-', false);
 position_3
 scatter3(alpha_3(1), position_3(1), position_3(1+Niter), 'green')
 
rotate on

%% Surface for the angle alpha_{n+1}

% Cobweb surface

epsilon = 0.001;

for side_id = 0:3

    [alpha_0 P_0] = meshgrid(0 + epsilon:0.05:pi-epsilon, side_id + epsilon:0.01:side_id+1-epsilon); 
    % Generate x = alpha_n and y = P_{n} data
    % Not good to include starting angles of exactly 0 and pi as these result
    % in discontinuties in our plot.

    alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

    for row=1:size(alpha_0,1) % over the rows

        for col=1:size(alpha_0,2) % over the columns
            [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 4, ...
                                                                                        'o', false); % no plot
                                                                                    
            ith_iterate = 2;   %depends on the N we have chosen , x_1, x_2, x_3,...      
            % change this for any value from 2 to N+1, 
            % for period ith_iterate - 1
            
            alpha_1(row, col) = alpha_temp(ith_iterate);
            P_1(row, col) = pos_temp(ith_iterate);
        end
    end
    
    % Converting into vector as scatter3 wants vectors not matrices
    alpha_0_vec = alpha_0(:); P_0_vec = P_0(:); alpha_1_vec = alpha_1(:);
    
    S = repmat([20],numel(alpha_0_vec),1);
    C = repmat([3],numel(alpha_0_vec),1);
    s = S(:);
    c = C(:);

    scatter3(alpha_0_vec, P_0_vec, alpha_1_vec, c, s)
    %surf(alpha_0, P_0, P_1, 'Linestyle', 'None') %Plot the surface
    %surf(alpha_0, P_0, P_1) %Plot the surface
    xlabel('\alpha_n')
    ylabel('P_n')
    zlabel(sprintf('\alpha_{n+%d}', ith_iterate-1))
    hold on
end

set(gca,'XLim',[0 pi], 'YLim', [0 4], 'ZLim',[0 pi])
xticks([0 pi/2 pi])
xticklabels({'0','\pi/2', '\pi'})
zticks([0 pi/2 pi])
zticklabels({'0','\pi/2', '\pi'})
title(sprintf('Period: %d', ith_iterate-1))

%[alpha_0 P_0] = meshgrid(0 + epsilon:0.05:pi-epsilon, 0 + epsilon:0.01:4+1-epsilon); 
%P_1 = P_0;
%plane = surf(alpha_0, P_0, P_1, 'FaceAlpha', 1, 'FaceColor', 'r')

[x y] = meshgrid(0:0.01:pi, 0:0.1:4); % Generate x = alpha_n and y = P_{n} data
z = x;     % z = P_{n+1}

s = surf(x,y,z,'FaceAlpha',0.5, 'FaceColor', 'r') %Plot the surface
s.EdgeColor = 'none';

%set(plane,'facealpha',0.2)

rotate on

%% It should be sufficient to just show the plot for side 0
% since all the sides on the square are the same

% However we also provide this plot for all 3

% Drawing a plane

% Cobweb surface

epsilon = 0.001;

for side_id = 0:3

    [alpha_0 P_0] = meshgrid(0 + epsilon:0.05:pi-epsilon, side_id + epsilon:0.01:side_id+1-epsilon); 
    % Generate x = alpha_n and y = P_{n} data
    % Not good to include starting angles of exactly 0 and pi as these result
    % in discontinuties in our plot.

    alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

    for row=1:size(alpha_0,1) % over the rows

        for col=1:size(alpha_0,2) % over the columns
            [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 5, ...
                                                                                        'o', false); % no plot
                                                                                    
            ith_iterate = 4;   %depends on the N we have chosen , x_1, x_2, x_3,...      
            % change this for any value from 2 to N+1, 
            % for a period  of ith_iterate - 1
            
            alpha_1(row, col) = alpha_temp(ith_iterate);
            P_1(row, col) = pos_temp(ith_iterate);
        end
    end
    
    % Converting into vector as scatter3 wants vectors not matrices
    alpha_0_vec = alpha_0(:); P_0_vec = P_0(:); P_1_vec = P_1(:);
    
    S = repmat([20],numel(alpha_0_vec),1);
    C = repmat([3],numel(alpha_0_vec),1);
    s = S(:);
    c = C(:);

    scatter3(alpha_0_vec, P_0_vec, P_1_vec, c, s)
    %surf(alpha_0, P_0, P_1, 'Linestyle', 'None') %Plot the surface
    %surf(alpha_0, P_0, P_1) %Plot the surface
    xlabel('\alpha_n')
    ylabel('P_n')
    zlabel(sprintf('P_{n+%d}', ith_iterate-1))
    hold on
end

set(gca,'XLim',[0 pi], 'YLim', [0 1], 'ZLim',[0 4]) % We just needed to change the gca
xticks([0 pi/2 pi])
xticklabels({'0','\pi/2', '\pi'})
title(sprintf('Period: %d', ith_iterate-1))

%[alpha_0 P_0] = meshgrid(0 + epsilon:0.05:pi-epsilon, 0 + epsilon:0.01:4+1-epsilon); 
%P_1 = P_0;
%plane = surf(alpha_0, P_0, P_1, 'FaceAlpha', 1, 'FaceColor', 'r')

[x y] = meshgrid(0:0.01:pi, 0:0.1:4); % Generate x = alpha_n and y = P_{n} data
z = y;     % z = P_{n+1}

s = surf(x,y,z,'FaceAlpha',0.5, 'FaceColor', 'r') %Plot the surface
s.EdgeColor = 'none';

%set(plane,'facealpha',0.2)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Investigating period 3 orbits

% A candidate (this initial angle and position are chosen for us to narrowly miss corner 2)
Niter = 3;

[side_3 alpha_3 position_3] = square_billiards_cobweb(pi/4+0.001, 0.001, Niter, 'o-', false);

scatter3(alpha_3(1), position_3(1), position_3(1+Niter), 'blue')

 % Another candidate (this initial angle and position are chosen for us to narrowly miss corner 2)
 [side_3 alpha_3 position_3] = square_billiards_cobweb(acot(1/2)+0.00001, 0.5, Niter, 'o-', false);
 position_3
 scatter3(alpha_3(1), position_3(1), position_3(1+Niter), 'blue')
 
% Another candidate on other surface
% (this initial angle and position are chosen for us to narrowly miss corner 3)
[side_3 alpha_3 position_3] = square_billiards_cobweb(pi - acot(1/2)+0.00001, 0.5, Niter, 'o-', false);
 position_3
 scatter3(alpha_3(1), position_3(1), position_3(1+Niter), 'green')
 
rotate on




%% Using two distinct color maps to provide stark constrast near the pi/2 area. This allows us 
% to exactly tell the difference between the pi/2-alpha_0 branch and the
% pi/2+alpha_0 branch.

% Our surface

epsilon = 0.001;

for side_id = 0:3

    [alpha_0 P_0] = meshgrid(0 + epsilon:0.01:pi-epsilon, side_id + epsilon:0.05:side_id+1-epsilon); 
    % Generate x = alpha_n and y = P_{n} data
    % Not good to include starting angles of exactly 0 and pi as these result
    % in discontinuties in our plot.

    alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

    for row=1:size(alpha_0,1) % over the rows

        for col=1:size(alpha_0,2) % over the columns
            [side_temp alpha_temp pos_temp] = square_billiards_cobweb(alpha_0(row,col), P_0(row,col), 5, ...
                                                                                        'o', false); % no plot
                                                                                    
            ith_iterate = 4;   %depends on the N we have chosen , x_1, x_2, x_3,...      
            % change this for any value from 2 to N+1, 
            % for a period  of ith_iterate - 1
            
            alpha_1(row, col) = alpha_temp(ith_iterate);
            P_1(row, col) = pos_temp(ith_iterate);
        end
    end
    
    % Converting into vector as scatter3 wants vectors not matrices
    alpha_0_vec = alpha_0(:); P_0_vec = P_0(:); P_1_vec = P_1(:);
    
    S = repmat([20],numel(alpha_0_vec),1);
    C = repmat([3],numel(alpha_0_vec),1);
    s = S(:);
    c = C(:);
    
    %surf(alpha_0, P_0, P_1, alpha_1, 'Linestyle', 'None'), hold on % No mesh on surface
    % color depends on alpha_1
    surface_no_wall_discont(alpha_0, P_0, P_1, alpha_1)
    
    cmap = [ cool(64); winter(64) ];
    colormap(gca, cmap);
    clims = [0 pi];  % minimum and maximum alpha values
    caxis(clims);
    
    c = colorbar('Limits', [0 pi], 'Ticks', [0 pi], 'TickLabels', {'0', '\pi'}, 'FontSize', 10)
    c.Label.String = '\alpha_{n+3}'; c.Label.FontSize = 15;
    c.Label.Rotation = 90; % to rotate the label
    
   
    xlabel('\alpha_n')
    ylabel('P_n')
    zlabel(sprintf('P_{n+%d}', ith_iterate-1))
    hold on
end

set(gca,'XLim',[0 pi], 'YLim', [0 1], 'ZLim',[0 4])
xticks([0 pi/2 pi])
xticklabels({'0','\pi/2', '\pi'})
title(sprintf('Period: %d', ith_iterate-1))


% Our cobweb plane, P_n = P_n+1
[x y] = meshgrid(0:0.01:pi, 0:0.1:4); % Generate x = alpha_n and y = P_{n} data
z = y;     % z = P_{n+1}

s = surf(x,y,z,'FaceAlpha', 0.5, 'FaceColor', 'r') %Plot the surface
s.EdgeColor = 'none';

%set(plane,'facealpha',0.2)

rotate on




%% Testing orbit diagrams

N = 1000;

plot_N = 800;     % we only want to plot long term behaviour, from 0 to 300 is effectively the burn in period

P_0 = 0.5; % Initial condition/position
epsilon = 0.001;

for alpha = [0 + epsilon:0.1:pi-epsilon]     % range of initial angles
    
    [side_1 alpha_1 position_1] = square_billiards_cobweb(alpha, P_0, N, 'o', false);
    plot(alpha, position_1(plot_N:N), 'ro', 'MarkerSize', 3), hold on
end

title(sprintf('Orbit diagram for square billiard with constant initial position %0.2f', P_0))
xlabel('alpha - intial angle')
ylabel(sprintf('P_{n} for  %d<n<%d', plot_N, N))
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0', 'pi/4', '\pi/2', '3pi/4', '\pi'})
set(gca,'XLim',[0 pi])

% None of the ones above are periodic (or at least a small period), however one may investigate the one
% near 0.001+26*0.1 further as there might be a periodic orbit nearby

% Notice the lack of symmetry, we would expect that the behaviour for alpha
% and pi-alpha to be identical.

%% Zooming into the behaviour for alpha = 2.6

N = 5000;   %increasing the number of iterates for more accurate long term behaviour

plot_N = 4700;

P_0 = 0.5; % Initial condition/position
epsilon = 0.001;

for alpha = [2.59 + epsilon:0.0003:2.61-epsilon]     % range of initial angles
    
    [side_1 alpha_1 position_1] = square_billiards_cobweb(alpha, P_0, N, 'o', false);
    plot(alpha, position_1(plot_N:N), 'ro', 'MarkerSize', 3), hold on
end

title(sprintf('Orbit diagram for square billiard with constant initial position %0.2f', P_0))
xlabel('alpha - intial angle')
ylabel(sprintf('P_{n} for  %d<n<%d', plot_N, N))
%xticks([0 pi/4 pi/2 3*pi/4 pi])
%xticklabels({'0', 'pi/4', '\pi/2', '3pi/4', '\pi'})
%set(gca,'XLim',[0 pi])

pause     % press any key to show the actual period 16 orbit that lies within this range
[side_26 alpha_26 position_26] = square_billiards_cobweb(pi-atan(3/5), P_0, N, 'o', false);
plot(pi-atan(3/5), position_26(plot_N:N), 'bo', 'MarkerSize', 3), hold on

position_26(4968:5000)    %checking to see if it is indeed periodic (period 16)
    
%% Investigating the strange seemingly also period 16 orbit. This is the 35th vertical line in the plot
% corresponds to the angle of 2.6015

% NOTE: the two only possible combinations of period 16 are 3/5 and 5/3 as the slopes
% (simplfied fraction)

[side_test alpha_test position_test] = square_billiards_cobweb(2.6015, P_0, N, 'o', false);

position_test(4968:5000) %checking to see if it is indeed periodic , it is close to period 16 but it is slightly off

%% As above but Orbit diagram for alpha, angle as y axis

for alpha = [0 + epsilon:0.1:pi-epsilon]     % range of initial angles
    [side_1 alpha_1 position_1] = square_billiards_cobweb(alpha, P_0, N, 'o', false);
    plot(alpha, alpha_1(plot_N:N), 'ro', 'MarkerSize', 3), hold on
end


title(sprintf('Orbit diagram for square billiard with constant initial position %0.2f', P_0))
xlabel('alpha - intial angle')
ylabel(sprintf('alpha_{n} for  %d<n<%d', plot_N, N))
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0', 'pi/4', '\pi/2', '3pi/4', '\pi'})
set(gca,'XLim',[0 pi])

%% Testing some rational slopes as angles

rational_slopes = 1./[1, 3:20]; 
periodic_alpha = atan(rational_slopes); % range of initial angles

for i = 1:size(periodic_alpha, 2)   %iterate over the given angles
    alpha = periodic_alpha(i);
    [side_1 alpha_1 position_1] = square_billiards_cobweb(alpha, P_0, N, 'o', false);
    plot(position_1(plot_N:N), alpha, 'bo', 'MarkerSize', 3), hold on; pbaspect([2 1 1])
end

title(sprintf('Orbit diagram for square billiard with constant initial position %0.2f', P_0))
ylabel('\alpha - intial angle')
xlabel(sprintf('P_{n} for  %d<n<%d', plot_N, N))
xticks([1 2 3 4])
xticklabels({'1' '2' '3' '4'})
yticks([0 pi/4 pi/2 3*pi/4 pi])
yticklabels({'0', 'pi/4', '\pi/2', '3pi/4', '\pi'})
set(gca, 'XLim',[1 4], 'YLim',[0 pi/2])

%% Orbit diagram surface by incorportating the angle alpha as well

rational_slopes = 1./[1:10];
periodic_alpha = atan(rational_slopes); % range of initial angles

for i = 1:size(periodic_alpha, 2)   %iterate over the given angles
    alpha = periodic_alpha(i);
    [side_1 alpha_1 position_1] = square_billiards_cobweb(alpha, P_0, N, 'o', false);
    plot(alpha, alpha_1(plot_N:N), 'ro', 'MarkerSize', 3), hold on
end

title(sprintf('Orbit diagram for square billiard with constant initial position %0.2f', P_0))
xlabel('alpha - intial angle')
ylabel(sprintf('alpha_{n} for  %d<n<%d', plot_N, N))
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0', 'pi/4', '\pi/2', '3pi/4', '\pi'})
set(gca,'XLim',[0 pi])

%% Zooming into the behaviour for alpha near pi/2

N = 5000;   %increasing the number of iterates for more accurate long term behaviour

plot_N = 4700;

P_0 = 0.5; % Initial condition/position
epsilon = 0.001;

for alpha = [(pi/2 - 0.05):0.001:(pi/2 + 0.05)]     % range of initial angles
    
    [side_1 alpha_1 position_1] = square_billiards_cobweb(alpha, P_0, N, 'o', false);
    plot(alpha, position_1(plot_N:N), 'ro', 'MarkerSize', 3), hold on
end

title(sprintf('Orbit diagram for square billiard with constant initial position %0.2f', P_0))
xlabel('alpha - intial angle')
ylabel(sprintf('P_{n} for  %d<n<%d', plot_N, N))
xticks([pi/2-0.05 pi/2 pi/2+0.05])
xticklabels({'pi/2-0.05', '\pi/2', 'pi/2+0.05'})
set(gca,'XLim',[pi/2-0.05 pi/2+0.05])

pause     % press any key to show the period 2 orbit of angle pi/2 in blue
[side_2 alpha_2 position_2] = square_billiards_cobweb(pi/2, P_0, N, 'o', false);
plot(pi/2, position_2(plot_N:N), 'bo', 'MarkerSize', 3), hold on

% It is expected that for alpha approx pi/2, the billiard tends to
% oscillate between sides 0 and side 2 (indicated by the two thick bands).
% Near pi/2 some strange behaviour ends up occuring as the orbits begin to
% not cover the entire sides (indicated by the white space) before finally
% approaching pi/2.

%% Zooming into the behaviour for alpha near pi/4

N = 5000;   %increasing the number of iterates for more accurate long term behaviour

plot_N = 4700;

P_0 = 0.5; % Initial condition/position
epsilon = 0.001;

for alpha = [(pi/4 - 0.05):0.001:(pi/4 + 0.05)]     % range of initial angles
    
    [side_1 alpha_1 position_1] = square_billiards_cobweb(alpha, P_0, N, 'o', false);
    plot(alpha, position_1(plot_N:N), 'ro', 'MarkerSize', 3), hold on
end

title(sprintf('Orbit diagram for square billiard with constant initial position %0.2f', P_0))
xlabel('alpha - intial angle')
ylabel(sprintf('P_{n} for  %d<n<%d', plot_N, N))
xticks([pi/4-0.05 pi/4 pi/4+0.05])
xticklabels({'pi/4-0.05', '\pi/4', 'pi/4+0.05'})
set(gca,'XLim',[pi/4-0.05 pi/4+0.05])

pause     % press any key to show the period 2 orbit of angle pi/4 in blue
[side_4 alpha_4 position_4] = square_billiards_cobweb(pi/4, P_0, N, 'o', false);
plot(pi/4, position_4(plot_N:N), 'bo', 'MarkerSize', 3), hold on


%% Orbit diagrams for different initial position

N = 1000;

plot_N = 800;     % we only want to plot long term behaviour, from 0 to 300 is effectively the burn in period

P_0 = 0.27; % Initial condition/position
epsilon = 0.001;

for alpha = [0 + epsilon:0.1:pi-epsilon]     % range of initial angles
    
    [side_1 alpha_1 position_1] = square_billiards_cobweb(alpha, P_0, N, 'o', false);
    plot(alpha, position_1(plot_N:N), 'ro', 'MarkerSize', 3), hold on
end

title(sprintf('Orbit diagram for square billiard with constant initial position %0.2f', P_0))
xlabel('alpha - intial angle')
ylabel(sprintf('P_{n} for  %d<n<%d', plot_N, N))
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0', 'pi/4', '\pi/2', '3pi/4', '\pi'})
set(gca,'XLim',[0 pi])

% Notice how similar it is to the initial condition of 0.5. This is
% because the initial position has no effect on the periodicity, only the
% angle.



%%
% SQUARE BILLIARDS PHASE PORTRAITS
epsilon = 1e-3;
init_pos = [0+epsilon:0.05:1-epsilon];

% Period 2 orbits
for i=1:length(init_pos)
    [side alpha position] = square_billiards_cobweb(pi/2, init_pos(i), 2, 'o-', false);
    plot(position, alpha, 'bo', 'MarkerFaceColor', 'b','MarkerSize', 2); hold on;
end

% Period 4 orbits
for i=1:length(init_pos)
    [side alpha position] = square_billiards_cobweb(pi/4, init_pos(i), 4, 'o-', false);
    plot(position, alpha, 'go', 'MarkerFaceColor', 'g','MarkerSize', 2); hold on;
end

% Period 6 orbits
for i=1:length(init_pos)
    [side alpha position] = square_billiards_cobweb(atan(2/1), init_pos(i), 6, 'o-', false);
    plot(position, alpha, 'mo', 'MarkerFaceColor', 'g','MarkerSize', 2); hold on;
end

% Period 8 orbits
for i=1:length(init_pos)
    [side alpha position] = square_billiards_cobweb(atan(3/1), init_pos(i), 8, 'o-', false);
    plot(position, alpha, 'yo', 'MarkerFaceColor', 'g','MarkerSize', 2); hold on;
end

% Chaotic orbit 1
[side alpha position] = square_billiards_cobweb(pi/4+0.4, 0.3, N, 'o-', false);
plot(position, alpha, 'ro','MarkerSize', 2); hold on

% Chaotic orbit 2
[side alpha position] = square_billiards_cobweb(0.3, 0.3, N, 'o-', false);
plot(position, alpha, 'co','MarkerSize', 2); hold on

% Trajectories that hit the corners/vertex
% From side 1 to vertex 0. alpha_n = acot(1-P_n) + pi;
side = 1; P_n = [side+epsilon:0.05:side+1-epsilon]; alpha = acot(1-P_n) + pi; % to make sure we have positive angles
plot(P_n, alpha, 'k'); hold on

% From side 2 to vertex 0. alpha_n = acot(3-P_n);
side = 2; P_n = [side+epsilon:0.05:side+1-epsilon]; alpha = acot(3-P_n);
plot(P_n, alpha, 'k'); hold on

% From side 2 to vertex 1. alpha_n = acot(2-P_n) + pi;
side = 2; P_n = [side+epsilon:0.05:side+1-epsilon]; alpha = acot(2-P_n)+pi; % to make sure we have positive angles
plot(P_n, alpha, 'k'); hold on

% From side 3 to vertex 1. alpha_n = acot(4-P_n);
side = 3; P_n = [side+epsilon:0.05:side+1-epsilon]; alpha = acot(4-P_n);
plot(P_n, alpha, 'k'); hold on

% From side 0 to vertex 2. alpha_n = acot(1-P_n);
side = 0; P_n = [side+epsilon:0.05:side+1-epsilon]; alpha = acot(1-P_n);
plot(P_n, alpha, 'k'); hold on

% From side 3 to vertex 2. alpha_n = acot(3-P_n)+pi;
side = 3; P_n = [side+epsilon:0.05:side+1-epsilon]; alpha = acot(3-P_n)+pi;
plot(P_n, alpha, 'k'); hold on

% From side 0 to vertex 3. alpha_n = acot(-P_n)+pi;
side = 0; P_n = [side+epsilon:0.05:side+1-epsilon]; alpha = acot(0-P_n)+pi;
plot(P_n, alpha, 'k'); hold on

% From side 1 to vertex 3. alpha_n = acot(2-P_n);
side = 1; P_n = [side+epsilon:0.05:side+1-epsilon]; alpha = acot(2-P_n);
plot(P_n, alpha, 'k'); hold on

title('Phase portrait for the square billiard')
set(gca,'XLim',[0 4], 'YLim',[0 pi])
xlabel('P - position')
xticks([0 1 2 3 4])
xticklabels({'0', '1', '2', '3', '4'})
ylabel('\alpha - angle')
yticks([0 pi/4 pi/2 3*pi/4 pi])
yticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'})


%% PRE-IMAGES OF ORBITS WHICH END UP IN THE CORNER/VERTEX

% We change the angle to be pi-alpha, this forces our trajectory to simply
% go in the opposite direction (reversing time).

N = 1; % how many iterations we go back

cmap = jet(8);

% From side 1 to vertex 0. alpha_n = acot(1-P_n) + pi;
side = 1; P_n = [side+epsilon:0.01:side+1-epsilon]; alpha = acot(1-P_n) + pi; % to make sure we have positive angles

for i=1:length(P_n)
    [pre_side pre_alpha pre_position] = square_billiards_cobweb(pi-alpha(i), P_n(i), N, 'o-', false); % note the pi-alpha
    plot(pre_position, pre_alpha, 'o','MarkerSize', 2, 'color', cmap(1,:), 'MarkerFaceColor', cmap(8,:)); hold on
end

% From side 2 to vertex 0. alpha_n = acot(3-P_n);
side = 2; P_n = [side+epsilon:0.01:side+1-epsilon]; alpha = acot(3-P_n);

for i=1:length(P_n)
    [pre_side pre_alpha pre_position] = square_billiards_cobweb(pi-alpha(i), P_n(i), N, 'o-', false);
    plot(pre_position, pre_alpha, 'o','MarkerSize', 2, 'color', cmap(2,:), 'MarkerFaceColor', cmap(8,:)); hold on
end

% From side 2 to vertex 1. alpha_n = acot(2-P_n) + pi;
side = 2; P_n = [side+epsilon:0.01:side+1-epsilon]; alpha = acot(2-P_n)+pi; % to make sure we have positive angles

for i=1:length(P_n)
    [pre_side pre_alpha pre_position] = square_billiards_cobweb(pi-alpha(i), P_n(i), N, 'o-', false);
    plot(pre_position, pre_alpha, 'o','MarkerSize', 2, 'color', cmap(3,:), 'MarkerFaceColor', cmap(8,:)); hold on
end

% From side 3 to vertex 1. alpha_n = acot(4-P_n);
side = 3; P_n = [side+epsilon:0.01:side+1-epsilon]; alpha = acot(4-P_n);

for i=1:length(P_n)
    [pre_side pre_alpha pre_position] = square_billiards_cobweb(pi-alpha(i), P_n(i), N, 'o-', false);
    plot(pre_position, pre_alpha, 'o','MarkerSize', 2, 'color', cmap(4,:), 'MarkerFaceColor', cmap(8,:)); hold on
end

% From side 0 to vertex 2. alpha_n = acot(1-P_n);
side = 0; P_n = [side+epsilon:0.01:side+1-epsilon]; alpha = acot(1-P_n);
for i=1:length(P_n)
    [pre_side pre_alpha pre_position] = square_billiards_cobweb(pi-alpha(i), P_n(i), N, 'o-', false);
    plot(pre_position, pre_alpha, 'o','MarkerSize', 2, 'color', cmap(5,:), 'MarkerFaceColor', cmap(8,:)); hold on
end

% From side 3 to vertex 2. alpha_n = acot(3-P_n)+pi;
side = 3; P_n = [side+epsilon:0.01:side+1-epsilon]; alpha = acot(3-P_n)+pi;
for i=1:length(P_n)
    [pre_side pre_alpha pre_position] = square_billiards_cobweb(pi-alpha(i), P_n(i), N, 'o-', false);
    plot(pre_position, pre_alpha, 'o','MarkerSize', 2, 'color', cmap(6,:), 'MarkerFaceColor', cmap(8,:)); hold on
end

% From side 0 to vertex 3. alpha_n = acot(-P_n)+pi;
side = 0; P_n = [side+epsilon:0.01:side+1-epsilon]; alpha = acot(0-P_n)+pi;
for i=1:length(P_n)
    [pre_side pre_alpha pre_position] = square_billiards_cobweb(pi-alpha(i), P_n(i), N, 'o-', false);
    plot(pre_position, pre_alpha, 'o','MarkerSize', 2, 'color', cmap(7,:), 'MarkerFaceColor', cmap(8,:)); hold on
end

% From side 1 to vertex 3. alpha_n = acot(2-P_n);
side = 1; P_n = [side+epsilon:0.01:side+1-epsilon]; alpha = acot(2-P_n);
for i=1:length(P_n)
    [pre_side pre_alpha pre_position] = square_billiards_cobweb(pi-alpha(i), P_n(i), N, 'o-', false);
    plot(pre_position, pre_alpha, 'o','MarkerSize', 2, 'color', cmap(8,:), 'MarkerFaceColor', cmap(8,:)); hold on
end

title('Pre-images for the square billiard')
set(gca,'XLim',[0 4], 'YLim',[0 pi])
xlabel('P - position')
xticks([0 1 2 3 4])
xticklabels({'0', '1', '2', '3', '4'})
ylabel('\alpha - angle')
yticks([0 pi/4 pi/2 3*pi/4 pi])
yticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'})


%% Combining pre-images and phase portrait

epsilon = 1e-3;

% Trajectories that hit the corners/vertex
% From side 1 to vertex 0. alpha_n = acot(1-P_n) + pi;
side = 1; P_n = [side+epsilon:0.05:side+1-epsilon]; alpha = acot(1-P_n) + pi; % to make sure we have positive angles
plot(P_n, alpha, 'k'); hold on

% From side 2 to vertex 0. alpha_n = acot(3-P_n);
side = 2; P_n = [side+epsilon:0.05:side+1-epsilon]; alpha = acot(3-P_n);
plot(P_n, alpha, 'k'); hold on

% From side 2 to vertex 1. alpha_n = acot(2-P_n) + pi;
side = 2; P_n = [side+epsilon:0.05:side+1-epsilon]; alpha = acot(2-P_n)+pi; % to make sure we have positive angles
plot(P_n, alpha, 'k'); hold on

% From side 3 to vertex 1. alpha_n = acot(4-P_n);
side = 3; P_n = [side+epsilon:0.05:side+1-epsilon]; alpha = acot(4-P_n);
plot(P_n, alpha, 'k'); hold on

% From side 0 to vertex 2. alpha_n = acot(1-P_n);
side = 0; P_n = [side+epsilon:0.05:side+1-epsilon]; alpha = acot(1-P_n);
plot(P_n, alpha, 'k'); hold on

% From side 3 to vertex 2. alpha_n = acot(3-P_n)+pi;
side = 3; P_n = [side+epsilon:0.05:side+1-epsilon]; alpha = acot(3-P_n)+pi;
plot(P_n, alpha, 'k'); hold on

% From side 0 to vertex 3. alpha_n = acot(-P_n)+pi;
side = 0; P_n = [side+epsilon:0.05:side+1-epsilon]; alpha = acot(0-P_n)+pi;
plot(P_n, alpha, 'k'); hold on

% From side 1 to vertex 3. alpha_n = acot(2-P_n);
side = 1; P_n = [side+epsilon:0.05:side+1-epsilon]; alpha = acot(2-P_n);
plot(P_n, alpha, 'k'); hold on

% PRE-IMAGES are with black dots
N_reverse = 6; % number of iterations backwards

% From side 1 to vertex 0. alpha_n = acot(1-P_n) + pi;
side = 1; P_n = [side+epsilon:0.01:side+1-epsilon]; alpha = acot(1-P_n) + pi; % to make sure we have positive angles

for i=1:length(P_n)
    [pre_side pre_alpha pre_position] = square_billiards_cobweb(pi-alpha(i), P_n(i), N_reverse, 'o-', false); % note the pi-alpha
    plot(pre_position, pre_alpha, 'ko','MarkerSize', 1, 'MarkerFaceColor', 'k'); hold on
end

% From side 2 to vertex 0. alpha_n = acot(3-P_n);
side = 2; P_n = [side+epsilon:0.01:side+1-epsilon]; alpha = acot(3-P_n);

for i=1:length(P_n)
    [pre_side pre_alpha pre_position] = square_billiards_cobweb(pi-alpha(i), P_n(i), N_reverse, 'o-', false);
    plot(pre_position, pre_alpha, 'ko','MarkerSize', 1, 'MarkerFaceColor', 'k'); hold on
end

% From side 2 to vertex 1. alpha_n = acot(2-P_n) + pi;
side = 2; P_n = [side+epsilon:0.01:side+1-epsilon]; alpha = acot(2-P_n)+pi; % to make sure we have positive angles

for i=1:length(P_n)
    [pre_side pre_alpha pre_position] = square_billiards_cobweb(pi-alpha(i), P_n(i), N_reverse, 'o-', false);
    plot(pre_position, pre_alpha, 'ko','MarkerSize', 1, 'MarkerFaceColor', 'k'); hold on
end

% From side 3 to vertex 1. alpha_n = acot(4-P_n);
side = 3; P_n = [side+epsilon:0.01:side+1-epsilon]; alpha = acot(4-P_n);

for i=1:length(P_n)
    [pre_side pre_alpha pre_position] = square_billiards_cobweb(pi-alpha(i), P_n(i), N_reverse, 'o-', false);
    plot(pre_position, pre_alpha, 'ko','MarkerSize', 1, 'MarkerFaceColor', 'k'); hold on
end

% From side 0 to vertex 2. alpha_n = acot(1-P_n);
side = 0; P_n = [side+epsilon:0.01:side+1-epsilon]; alpha = acot(1-P_n);
for i=1:length(P_n)
    [pre_side pre_alpha pre_position] = square_billiards_cobweb(pi-alpha(i), P_n(i), N_reverse, 'o-', false);
    plot(pre_position, pre_alpha, 'ko','MarkerSize', 1, 'MarkerFaceColor', 'k'); hold on
end

% From side 3 to vertex 2. alpha_n = acot(3-P_n)+pi;
side = 3; P_n = [side+epsilon:0.01:side+1-epsilon]; alpha = acot(3-P_n)+pi;
for i=1:length(P_n)
    [pre_side pre_alpha pre_position] = square_billiards_cobweb(pi-alpha(i), P_n(i), N_reverse, 'o-', false);
    plot(pre_position, pre_alpha, 'ko','MarkerSize', 1, 'MarkerFaceColor', 'k'); hold on
end

% From side 0 to vertex 3. alpha_n = acot(-P_n)+pi;
side = 0; P_n = [side+epsilon:0.01:side+1-epsilon]; alpha = acot(0-P_n)+pi;
for i=1:length(P_n)
    [pre_side pre_alpha pre_position] = square_billiards_cobweb(pi-alpha(i), P_n(i), N_reverse, 'o-', false);
    plot(pre_position, pre_alpha, 'ko','MarkerSize', 1, 'MarkerFaceColor', 'k'); hold on
end

% From side 1 to vertex 3. alpha_n = acot(2-P_n);
side = 1; P_n = [side+epsilon:0.01:side+1-epsilon]; alpha = acot(2-P_n);
for i=1:length(P_n)
    [pre_side pre_alpha pre_position] = square_billiards_cobweb(pi-alpha(i), P_n(i), N_reverse, 'o-', false);
    plot(pre_position, pre_alpha, 'ko','MarkerSize', 1, 'MarkerFaceColor', 'k'); hold on
end


init_pos = [0+epsilon:0.05:1-epsilon];

% Period 2 orbits
for i=1:length(init_pos)
    [side alpha position] = square_billiards_cobweb(pi/2, init_pos(i), 2, 'o-', false);
    plot(position, alpha, 'bo', 'MarkerFaceColor', 'b','MarkerSize', 2); hold on;
end

% Period 4 orbits
for i=1:length(init_pos)
    [side alpha position] = square_billiards_cobweb(pi/4, init_pos(i), 4, 'o-', false);
    plot(position, alpha, 'go', 'MarkerFaceColor', 'g','MarkerSize', 2); hold on;
end

% Period 6 orbits
for i=1:length(init_pos)
    [side alpha position] = square_billiards_cobweb(atan(2/1), init_pos(i), 6, 'o-', false);
    plot(position, alpha, 'mo', 'MarkerFaceColor', 'g','MarkerSize', 2); hold on;
end

% Period 8 orbits
for i=1:length(init_pos)
    [side alpha position] = square_billiards_cobweb(atan(3/1), init_pos(i), 8, 'o-', false);
    plot(position, alpha, 'yo', 'MarkerFaceColor', 'g','MarkerSize', 2); hold on;
end

N = 1000;
% Chaotic orbit 1
[side alpha position] = square_billiards_cobweb(pi/4+0.4, 0.3, N, 'o-', false);
plot(position, alpha, 'ro','MarkerSize', 2); hold on

% Chaotic orbit 2
[side alpha position] = square_billiards_cobweb(0.3, 0.3, N, 'o-', false);
plot(position, alpha, 'co','MarkerSize', 2); hold on


title('Phase portrait for the square billiard')
set(gca,'XLim',[0 4], 'YLim',[0 pi])
xlabel('P - position')
xticks([0 1 2 3 4])
xticklabels({'0', '1', '2', '3', '4'})
ylabel('\alpha - angle')
yticks([0 pi/4 pi/2 3*pi/4 pi])
yticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'})


%% TESTING COLOR MAPS

epsilon = 1e-4;
init_angle = [0+epsilon:0.01:pi/2-epsilon];

cmap = jet(length(init_angle));

for i=1:length(init_angle)
    [side alpha position] = square_billiards_cobweb(init_angle(i), 0.3, N, 'o-', false);
    plot(position, alpha, 'o', 'color', cmap(i,:),'MarkerSize', 1, 'MarkerFaceColor', cmap(i,:)); hold on;
end

colormap(gca, cmap);
clims = [0 pi/2];  % minimum and maximum alpha values
caxis(clims);
c = colorbar('Limits', [0 pi/2], 'Ticks', [0 pi/2], 'TickLabels', {'0', '\pi/2'}, 'FontSize', 10)
c.Label.String = '\alpha_{0}'; c.Label.FontSize = 15;
c.Label.Rotation = 0; % to rotate the label
%pos = get(c,'Position');
%c.Label.Position = [pos(1)+1.5 pos(2)+pi/2]; % to change its position
    
title('Phase portrait for the square billiard')
set(gca,'XLim',[0 4], 'YLim',[0 pi])
xlabel('P - position')
xticks([0 1 2 3 4])
xticklabels({'0', '1', '2', '3', '4'})
ylabel('\alpha - angle')
yticks([0 pi/4 pi/2 3*pi/4 pi])
yticklabels({'0', 'pi/4', '\pi/2', '3\pi/4', '\pi'})


%% P6 Transcritical/Period doubling

% Period 6 orbits
init_pos_first_half = [0+epsilon:0.01:0.5-epsilon];
init_pos = [0+epsilon:0.01:4-epsilon];

%init_pos = [0.3]
for i=1:length(init_pos)
    [side alpha position] = square_billiards_cobweb(atan(2/1), init_pos(i), 6, 'o-', false);
    plot(init_pos(i).*ones(length(position)), position, 'mo', 'MarkerFaceColor', 'g','MarkerSize', 2); hold on; pbaspect([1 1 1])
end

%plot(init_pos, 3-2*init_pos, 'b')
%plot(init_pos, 3/2+init_pos, 'g') % Derived analytically

%
% figuring out equation for the lines
%plot(init_pos, 2+2*init_pos, 'k')
%plot(init_pos_first_half, 4-2*init_pos_first_half, 'k'); hold on;
%plot(init_pos_first_half, 1+2*init_pos_first_half, 'k'); hold on;
%plot(init_pos_first_half, 5/2+init_pos_first_half, 'g'); hold on;  
%plot(init_pos_first_half, 5/2-init_pos_first_half, 'g'); hold on;
%plot(init_pos_first_half, init_pos_first_half, 'b'); hold on;
%plot(init_pos_first_half, 1-init_pos_first_half, 'b')

set(gca,'XLim',[0 4], 'YLim',[0 4])
xlabel('P_0 - initial position')
xticks([0 1 2 3 4])
xticklabels({'0', '1', '2', '3', '4'})
ylabel('P_n - position')
yticks([0 1 2 3 4])
yticklabels({'0', '1', '2', '3', '4'})

%% TEST 

init_pos_first_half = [0+epsilon:0.01:0.5-epsilon];
init_pos = [0+epsilon:0.01:4-epsilon];

%init_pos = [0.3]
for i=1:length(init_pos)
    [side alpha position] = square_billiards_cobweb(atan(2/1), init_pos(i), 6, 'o-', false);
    plot3(init_pos(i).*ones(length(position)), alpha, position, 'mo', 'MarkerFaceColor', 'g','MarkerSize', 2); hold on;
    plot3(init_pos(i), alpha(1), position(1), 'bo','MarkerSize', 2)
end

set(gca,'XLim',[0 1], 'YLim',[0 4]);
xlabel('P_0 - initial position')
%xticks([0 1 2 3 4])
%xticklabels({'0', '1', '2', '3', '4'})

ylabel('\alpha_n - angle')

zlabel('P_n - position')
%yticks([0 pi/4 pi/2 3*pi/4 pi])
%yticklabels({'0', 'pi/4', '\pi/2', '3\pi/4', '\pi'})
rotate on;

%% P8 Transcritical/Period doubling

% Period 8 orbits
init_pos_first_half = [0+epsilon:0.01:0.5-epsilon];
init_pos = [0+epsilon:0.01:4-epsilon];

for i=1:length(init_pos)
    [side alpha position] = square_billiards_cobweb(atan(3/1)+0.0, init_pos(i), 50, 'o-', false);
    plot(init_pos(i).*ones(length(position)), position, 'go', 'MarkerFaceColor', 'g','MarkerSize', 2); hold on;
end

set(gca,'XLim',[0 4], 'YLim',[0 4])
xlabel('P_0 - initial position')
%xticks([0 1 2 3 4])
%xticklabels({'0', '1', '2', '3', '4'})
ylabel('P_n - position')
%yticks([0 pi/4 pi/2 3*pi/4 pi])
%yticklabels({'0', 'pi/4', '\pi/2', '3\pi/4', '\pi'})


%% P10 Transcritical/Period doubling

% Period 10 orbits
init_pos_first_half = [0+epsilon:0.01:0.5-epsilon];
init_pos = [0+epsilon:0.01:4-epsilon];

for i=1:length(init_pos)
    [side alpha position] = square_billiards_cobweb(atan(4)+0.0, init_pos(i), 15, 'o-', false);
    plot(init_pos(i).*ones(length(position)), position, 'bo', 'MarkerFaceColor', 'g','MarkerSize', 2); hold on;  pbaspect([1 1 1])
end

set(gca,'XLim',[0 1], 'YLim',[0 1])
xlabel('P_0 - initial position')
%xticks([0 1 2 3 4])
%xticklabels({'0', '1', '2', '3', '4'})
ylabel('P_n - position')
%yticks([0 1 2 3 4])
%yticklabels({'0', '1', '2', '3', '4'})



%% P14 Transcritical/Period doubling

% Period 14 orbits
init_pos = [0+epsilon:0.005:1-epsilon];

for i=1:length(init_pos)
    [side alpha position] = square_billiards_cobweb(atan(5/2)+0.0, init_pos(i), 15, 'o-', false);
    plot(init_pos(i).*ones(length(position)), position, 'bo', 'MarkerFaceColor', 'g','MarkerSize', 2); hold on;  pbaspect([1 1 1])
end

set(gca,'XLim',[0 1], 'YLim',[0 1])
xlabel('P_0 - initial position')
xticks([0 1])
xticklabels({'0', '1'})
ylabel('P_n - position')
yticks([0 1])
yticklabels({'0', '1'})












%%

% Testing chaotic
init_pos_first_half = [0+epsilon:0.01:0.5-epsilon];
init_pos = [0+epsilon:0.01:4-epsilon];

%init_pos = [0.3]
for i=1:length(init_pos)
    [side alpha position] = square_billiards_cobweb(atan(2/1)+0.01, init_pos(i), 50, 'o-', false);
    %[side alpha position] = square_billiards_cobweb(1.5, init_pos(i), 10, 'o-', false);
    plot(init_pos(i).*ones(length(position)), position, 'mo', 'MarkerFaceColor', 'g','MarkerSize', 2); hold on;
end

set(gca,'XLim',[0 1], 'YLim',[0 4])
xlabel('P_0 - initial position')
%xticks([0 1 2 3 4])
%xticklabels({'0', '1', '2', '3', '4'})
ylabel('P_n - position')
%yticks([0 pi/4 pi/2 3*pi/4 pi])
%yticklabels({'0', 'pi/4', '\pi/2', '3\pi/4', '\pi'})