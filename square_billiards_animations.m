clear all; close all; clc; 

% Animations of the billiard trajectory

% Our Difference Equation for Square Billiards:
% x(n+1) = f(x(n)). f is too complex to write here. Refer to written notes.

N = 100;   % number of iterations/maps (bounces)

%% Period 4 orbit with lines and dots. Trivial angle of pi/4 and position of 0.5. 
% Follow the mapping from each position (dot in image)

init_pos = 0.2; init_angle = pi/8;
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'o-', false);
P_anim = square_billiards_draw(side, alpha, position);

%% Period 6 orbit

init_pos = 0.05; init_angle = atan(2/1);
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'ro-', false);
P_anim = square_billiards_draw(side, alpha, position);

%% Period 8 orbit

init_pos =0.5; % 2/3-0.03; 
init_angle = atan(3/1);
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'ro-', false);
P_anim = square_billiards_draw(side, alpha, position);




%% TEMP
N=10;
init_pos = 0.4; init_angle = atan(2);
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'ro-', false);
P_anim = square_billiards_draw(side, alpha, position);



%% Chaotic testing
% DELETE THIS 

init_pos = 0.929; init_angle = 1.5;
%init_pos = 0; init_angle = atan(7/5);
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, 100, 'ro-', false);
P_anim = square_billiards_draw(side, alpha, position);


%% Chaotic billiard (non-rational slope)

init_pos = 0.5; init_angle = atan(0.3235645);
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'o-', false);
P_anim = square_billiards_draw(side, alpha, position);

%% Period 120 orbit

init_pos = 0.5; init_angle = atan(23/37);
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'o-', false);
P_anim = square_billiards_draw(side, alpha, position);




%% TESTING

init_pos = 0.1; init_angle = atan(4);
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'ro-', false);
P_anim = square_billiards_draw(side, alpha, position);


%% With animations
% Period 4 starting from P=0.5
init_pos = 0.5; init_angle = pi/4;
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'o-', false);
P_anim = square_billiards_draw(side, alpha, position, 'animate');

%% Period 4 starting from P=0.2
init_pos = 0.2; init_angle = pi/4;
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'o-', false);
P_anim = square_billiards_draw(side, alpha, position, 'animate');

%% Period 6 starting from P=0.8
init_pos = 0.8; init_angle = atan(2);
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'o-', false);
P_anim = square_billiards_draw(side, alpha, position, 'animate');

%% Chaotic orbit
init_pos = 0.8; init_angle = atan(2)+0.1;
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'o-', false);
P_anim = square_billiards_draw(side, alpha, position, 'animate');



%% Triangle for appendix, proof that period 3 orbits are impossible in triangle

% Triangle with vertices at (0, 0.2), (1, 0.6) and (0.6, 1).
% The three equations for the sides/lines are:
% y = 0.4x+0.2, y = 4/3x+0.2, y = -x + 1.6

% The corners of the triangle as red dots
plot([0 0.6 1], [0.2 1 0.6], 'ro', 'linewidth', 2), hold on

% Connecting the corners with thick black lines (the billiard trajcetory)
x = [0:0.1:1];
plot(x, 0.4*x+0.2,'k', 'linewidth', 2), hold on
plot(x, 4/3*x+0.2,'k', 'linewidth', 2), hold on
plot(x, -x+1.6,'k', 'linewidth', 2)

% Labelling the vertices P, Q and R
text(-0.09, 0.2, 'R', 'FontSize', 22)
text(0.6, 1.05, 'Q', 'FontSize', 22)
text(1.03, 0.6, 'P', 'FontSize', 22)

% Angle Arcs
% For R
ang([0 0.2], 0.15, [0.95 pi/2],'r-')
ang([0 0.2], 0.12, [0.4 -pi/2],'r-')
text(0.1, 0.1, '\beta', 'Color', 'r', 'FontSize', 22)
text(0.03, 0.42, '\beta', 'Color', 'r', 'FontSize', 22)

% For P
ang([1 0.6], 0.15, [pi/2 pi/2+0.8],'r-')
ang([1 0.6], 0.12, [pi + 0.35 3*pi/2],'r-')
text(0.9, 0.45, '\alpha', 'Color', 'r', 'FontSize', 22)
text(0.92, 0.81, '\alpha', 'Color', 'r', 'FontSize', 22)

% For Q
ang([0.6 1], 0.15, [0 -0.77],'r-')
ang([0.6 1], 0.12, [pi pi + 0.9],'r-')
text(0.75, 0.92, '\theta', 'Color', 'r', 'FontSize', 22)
text(0.43, 0.92, '\theta', 'Color', 'r', 'FontSize', 22)

set(gca,'XLim',[0 1], 'YLim', [0 1])

% Add axes labels
ax1 = gca;
set(gca, 'linewidth', 3, 'FontSize', 22);
pbaspect([1 1 1])
set(gcf,'color','w');  % background color to white instead of usual grey
ax1.XColor = 'b';  ax1.YColor = 'b'; % change colour of axis

xticks(ax1, [])
xticklabels(ax1, {''})
yticks(ax1, [0 1])
yticklabels(ax1, {'A', 'D'})

%ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top',... 
%  'YAxisLocation','right','Color','none','XColor','b','YColor','b', ...
%  'linewidth', 3); % to plot on other two axes as well
%set(ax2, 'FontSize', 22)

%xticks(ax2, [])
%xticklabels(ax2, {''})
%yticks(ax2, [0 1])
%yticklabels(ax2, {'B', 'C'})

text(1.03, 0, 'B', 'Color','b', 'FontSize', 22)
text(1.03, 1, 'C', 'Color','b', 'FontSize', 22)

%% Plotting an "unfolding". A orbit with intial angle atan(1/2)+0.1 and position 0.2.

N = 6; % First reflection
init_angle = atan(1/2)+0.1; init_pos = 0.2;
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'o-', false);
P_anim = square_billiards_draw(side, alpha, position); % Convert into our square animation
text(-0.2, 1.05, '(a)', 'FontSize', 22) % label for the dissertation

%% Reflecting the table plot ONCE
plot(P_anim(1:2,1), P_anim(1:2,2), 'k', 'linewidth', 2), hold on % solid black for first line
plot(P_anim(2:3,1), P_anim(2:3,2), '--k', 'linewidth', 2), hold on % dotted black for reflected line

% New line that is in the second square
plot([2;2] - P_anim(2:3,1), P_anim(2:3,2), 'k', 'linewidth', 2), hold on

% Now the axis is not the sides of the square. We need to manually each of the sides.
% First square
plot([0 1], [0 0], 'b', 'linewidth', 3) % note we just need the vertices
plot([1 1], [0 1], 'b', 'linewidth', 3)
plot([0 0], [0 1], 'b', 'linewidth', 3)
plot([0 1], [1 1], 'b', 'linewidth', 3)

% Second square
plot([1 2], [0 0], 'b', 'linewidth', 3) % note we just need the vertices
plot([2 2], [0 1], 'b', 'linewidth', 3)
plot([1 1], [0 1], 'b', 'linewidth', 3)
plot([1 2], [1 1], 'b', 'linewidth', 3)

pbaspect([2 1 1])  % aspect ratio
set(gca,'XLim',[0 2], 'YLim', [0 1])

% Add label for the first square
text(-0.1, -0.1, 'A', 'Color','b', 'FontSize', 22)
text(1, -0.1, 'B', 'Color','b', 'FontSize', 22)
text(1, 1.1, 'C', 'Color','b', 'FontSize', 22)
text(-0.1, 1.1, 'D', 'Color','b', 'FontSize', 22)

% Add label for the second square
text(2, -0.1, 'A', 'Color','b', 'FontSize', 22)
text(2, 1.1, 'D', 'Color','b', 'FontSize', 22)

% Controlling the plot
set(gca, 'linewidth', 3, 'FontSize', 22);
set(gcf,'color','w');  % background color to white instead of usual grey
set(gca,'xtick', [], 'ytick', [], 'visible','off')  % removes axes and ticks

text(-0.3, 1.3, '(a)', 'FontSize', 22) %for the dissertation

%% Reflecting the table plot TWICE
%plot(P_anim(1:2,1), P_anim(1:2,2), 'k', 'linewidth', 2), hold on % solid black for first line
plot(P_anim(2:3,1), P_anim(2:3,2), '--k', 'linewidth', 2), hold on % dotted black for first reflected line
plot([2; 2]-P_anim(3:4,1), P_anim(3:4,2), '--k', 'linewidth', 2), hold on % second reflected line
plot(P_anim(4:5,1)+[2;0], [2;2]-P_anim(4:5,2), '--k', 'linewidth', 2), hold on % third reflected line
plot(P_anim(5:6,1)+[2;2], [2;2]-P_anim(5:6,2), '--k', 'linewidth', 2), hold on % fourth reflected line
plot([4;4]-P_anim(6:7,1), [2;2]-P_anim(6:7,2), '--k', 'linewidth', 2), hold on % fifth reflected line

% Now the axis is not the sides of the square. We need to manually each of the sides.
% First square
plot([0 1], [0 0], 'b', 'linewidth', 3) % note we just need the vertices
plot([1 1], [0 1], 'b', 'linewidth', 3)
plot([0 0], [0 1], 'b', 'linewidth', 3)
plot([0 1], [1 1], 'b', 'linewidth', 3)

% Second square
plot([1 2], [0 0], 'b', 'linewidth', 3) % add [1 1] to all the x co-ords of the first square
plot([2 2], [0 1], 'b', 'linewidth', 3) % since horizontal reflection
plot([1 1], [0 1], 'b', 'linewidth', 3)
plot([1 2], [1 1], 'b', 'linewidth', 3)

% Third square
plot([1 2], [1 1], 'b', 'linewidth', 3) % add [1 1] to all the y co-ords of the first square
plot([2 2], [1 2], 'b', 'linewidth', 3) % since vertical reflection
plot([1 1], [1 2], 'b', 'linewidth', 3)
plot([1 2], [2 2], 'b', 'linewidth', 3)

% Fourth square
plot([2 3], [1 1], 'b', 'linewidth', 3) % horizontal reflection
plot([3 3], [1 2], 'b', 'linewidth', 3)
plot([2 2], [1 2], 'b', 'linewidth', 3)
plot([2 3], [2 2], 'b', 'linewidth', 3)

% Fifth square
plot([3 4], [1 1], 'b', 'linewidth', 3) % horizontal reflection
plot([4 4], [1 2], 'b', 'linewidth', 3)
plot([3 3], [1 2], 'b', 'linewidth', 3)
plot([3 4], [2 2], 'b', 'linewidth', 3)

% Sixth square
plot([3 4], [2 2], 'b', 'linewidth', 3) % vertical reflection
plot([4 4], [2 3], 'b', 'linewidth', 3)
plot([3 3], [2 3], 'b', 'linewidth', 3)
plot([3 4], [3 3], 'b', 'linewidth', 3)


% Straightline trajectory
x = [0:5]; y = tan(init_angle)*x - init_angle * init_pos;
plot(x, y, 'k', 'linewidth', 2)

pbaspect([4 3 1])  % aspect ratio
set(gca,'XLim',[0 4], 'YLim', [0 3])

% Add label for the first square
text(-0.1, -0.15, 'A', 'Color','b', 'FontSize', 22)
text(0.75, -0.15, 'B', 'Color','b', 'FontSize', 22)
text(0.75, 1.15, 'C', 'Color','b', 'FontSize', 22)
text(-0.1, 1.15, 'D', 'Color','b', 'FontSize', 22)

% Add label for the second square
text(2.1, -0.15, 'A', 'Color','b', 'FontSize', 22)
text(2.1, 0.82, 'D', 'Color','b', 'FontSize', 22)

% Add label for the third square
text(0.75, 2.15, 'B', 'Color','b', 'FontSize', 22)
text(2.1, 2.15, 'A', 'Color','b', 'FontSize', 22)

% Add label for the fourth square
text(2.75, 2.15, 'B', 'Color','b', 'FontSize', 22)
text(2.75, 0.82, 'C', 'Color','b', 'FontSize', 22)

% Add label for the fifth square
text(4.1, 0.82, 'D', 'Color','b', 'FontSize', 22)
text(4.1, 2.15, 'A', 'Color','b', 'FontSize', 22)

% Add label for the sixth square
text(4.1, 3.15, 'D', 'Color','b', 'FontSize', 22)
text(2.75, 3.15, 'C', 'Color','b', 'FontSize', 22)

% Controlling the plot
set(gca, 'linewidth', 3, 'FontSize', 22);
set(gcf,'color','w');  % background color to white instead of usual grey
set(gca,'xtick', [], 'ytick', [], 'visible','off')  % removes axes and ticks

text(-0.3, 3.2, '(b)', 'FontSize', 22) %for the dissertation



%% Torus plot (2x2 square grid)

N = 6; % First reflection
init_angle = atan(1/2)+0.1; init_pos = 0.2;
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'o-', false);
P_anim = square_billiards_draw(side, alpha, position); % Convert into our square animation
text(-0.2, 1.05, '(a)', 'FontSize', 22) % label for the dissertation

%% Reflecting the table plot ONCE
plot(P_anim(1:2,1), P_anim(1:2,2), 'w', 'linewidth', 2), hold on % do not show the original lines (white color)
plot(P_anim(2:3,1), P_anim(2:3,2), '--w', 'linewidth', 2), hold on % (white color)

% New line that is in the second square
plot([2;2] - P_anim(2:3,1), P_anim(2:3,2), 'w', 'linewidth', 2), hold on

% Now the axis is not the sides of the square. We need to manually each of the sides.
% Original square
plot([0 1], [0 0], 'b', 'linewidth', 3) % note we just need the vertices
plot([1 1], [0 1], 'b', 'linewidth', 3)
plot([0 0], [0 1], 'b', 'linewidth', 3)
plot([0 1], [1 1], 'b', 'linewidth', 3)

% Right square
plot([1 2], [0 0], 'b', 'linewidth', 3) % note we just need the vertices
plot([2 2], [0 1], 'b', 'linewidth', 3)
plot([1 1], [0 1], 'b', 'linewidth', 3)
plot([1 2], [1 1], 'b', 'linewidth', 3)

% Diagonal square
plot([1 2], [1 1], 'b', 'linewidth', 3) % add [1 1] to all the y co-ords of the first square
plot([2 2], [1 2], 'b', 'linewidth', 3) % since vertical reflection
plot([1 1], [1 2], 'b', 'linewidth', 3)
plot([1 2], [2 2], 'b', 'linewidth', 3)

% Top square
plot([0 1], [2 2], 'b', 'linewidth', 3) 
plot([0 0], [1 2], 'b', 'linewidth', 3) 

pbaspect([1 1 1])  % aspect ratio
set(gca,'XLim',[0 2], 'YLim', [0 2])

% Add label for the inital square
text(-0.1, -0.1, 'A', 'Color','b', 'FontSize', 22)
text(1.03, -0.1, 'B', 'Color','b', 'FontSize', 22)
text(1.03, 1.1, 'C', 'Color','b', 'FontSize', 22)
text(-0.1, 1.1, 'D', 'Color','b', 'FontSize', 22)

% Add label for the right square
text(2.03, -0.1, 'A', 'Color','b', 'FontSize', 22)
text(2.03, 1.1, 'D', 'Color','b', 'FontSize', 22)

% Add label for the top square
text(-0.1, 2.1, 'A', 'Color','b', 'FontSize', 22)
text(1.03, 2.1, 'B', 'Color','b', 'FontSize', 22)

%Label for diagonal square
text(2.03, 2.1, 'A', 'Color','b', 'FontSize', 22)

% Drawing the lines:
% y=-0.3+0.4x, y=0.5+0.4x, y=1.25+0.4x
x = [0 2]; 
plot(x, -0.3+0.4*x, 'k', 'linewidth', 2);
plot(x, 0.5+0.4*x, 'k', 'linewidth', 2)
plot(x, 1.25+0.4*x, 'k', 'linewidth', 2)

% Controlling the plot
set(gca, 'linewidth', 3, 'FontSize', 22);
set(gcf,'color','w');  % background color to white instead of usual grey
set(gca,'xtick', [], 'ytick', [], 'visible','off')  % removes axes and ticks

text(-0.3, 2.1, '(b)', 'FontSize', 22) %for the dissertation



%% MY ALTERNATIVE PROOF FOR PERIODICTY
% Adapted from the earlier unfolding plot

%% Plotting an "unfolding". A orbit with intial angle atan(1/2) and position 0.2.

N = 6; % First reflection
init_angle = atan(1/2); init_pos = 0.2;
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'o-', false);
P_anim = square_billiards_draw(side, alpha, position); % Convert into our square animation
text(-0.2, 1.05, '(a)', 'FontSize', 22) % label for the dissertation


%% The unfolded plot

% Now the axis is not the sides of the square. We need to manually each of the sides.
% First square
plot([0 1], [0 0], 'b', 'linewidth', 3); hold on % note we just need the vertices
plot([1 1], [0 1], 'b', 'linewidth', 3)
plot([0 0], [0 1], 'b', 'linewidth', 3)
plot([0 1], [1 1], 'b', 'linewidth', 3)

% Second square (shift all x-cords by 1)
plot([1 2], [0 0], 'b', 'linewidth', 3) % add [1 1] to all the x co-ords of the first square
plot([2 2], [0 1], 'b', 'linewidth', 3) % since horizontal reflection
plot([1 1], [0 1], 'b', 'linewidth', 3)
plot([1 2], [1 1], 'b', 'linewidth', 3)

% Third square (shift all x-cords by 1)
plot([2 3], [0 0], 'b', 'linewidth', 3) % horizontal reflection
plot([3 3], [0 1], 'b', 'linewidth', 3)

% Fourth square (shift all y-cords by 1)
plot([2 3], [1 1], 'b', 'linewidth', 3) % vertical reflection
plot([3 3], [1 2], 'b', 'linewidth', 3)
plot([2 2], [1 2], 'b', 'linewidth', 3)
plot([2 3], [2 2], 'b', 'linewidth', 3)

% Fifth square (shift all x-cords by 1)
plot([3 4], [1 1], 'b', 'linewidth', 3) % horizontal reflection
plot([4 4], [1 2], 'b', 'linewidth', 3)
plot([3 3], [1 2], 'b', 'linewidth', 3)
plot([3 4], [2 2], 'b', 'linewidth', 3)

% Sixth square (shift all x-cords by 1)
plot([4 5], [2 2], 'b', 'linewidth', 3) % horizontal reflection
plot([4 5], [1 1], 'b', 'linewidth', 3)
plot([4 4], [1 2], 'b', 'linewidth', 3)
plot([5 5], [1 2], 'b', 'linewidth', 3)

% Seventh square (shift all y-cords by 1)
plot([4 5], [2 2], 'b', 'linewidth', 3) % vertical reflection
plot([5 5], [2 3], 'b', 'linewidth', 3)
plot([4 4], [2 3], 'b', 'linewidth', 3)
plot([4 5], [3 3], 'b', 'linewidth', 3)

% Eight square (shift all x-cords by 1)
plot([5 6], [2 2], 'b', 'linewidth', 3) % vertical reflection
plot([6 6], [2 3], 'b', 'linewidth', 3)
plot([5 5], [2 3], 'b', 'linewidth', 3)
plot([5 6], [3 3], 'b', 'linewidth', 3)

% Ninth square (shift all x-cords by 1)
%plot([6 7], [2 2], 'b', 'linewidth', 3) % vertical reflection
%plot([7 7], [2 3], 'b', 'linewidth', 3)
%plot([6 6], [2 3], 'b', 'linewidth', 3)
%plot([6 7], [3 3], 'b', 'linewidth', 3)

% Tenth square (shift all y-cords by 1)
%plot([6 7], [3 3], 'b', 'linewidth', 3) % vertical reflection
%plot([7 7], [3 4], 'b', 'linewidth', 3)
%plot([6 6], [3 4], 'b', 'linewidth', 3)
%plot([6 7], [4 4], 'b', 'linewidth', 3)

% Straightline trajectory
x = [0:7]; y = tan(init_angle)*x - init_angle * init_pos;
plot(x, y, 'k', 'linewidth', 2)

pbaspect([6 4 1])  % aspect ratio
set(gca,'XLim',[0 6], 'YLim', [0 4])

% Add label for the first square
text(-0.1, -0.15, 'A', 'Color','b', 'FontSize', 22)
text(0.75, -0.15, 'B', 'Color','b', 'FontSize', 22)
text(0.75, 1.15, 'C', 'Color','b', 'FontSize', 22)
text(-0.1, 1.15, 'D', 'Color','b', 'FontSize', 22)

% Add label for the second square
text(1.75, -0.15, 'A', 'Color','b', 'FontSize', 22)
text(1.75, 1.15, 'D', 'Color','b', 'FontSize', 22)

% Add label for the third square
text(3.05, 0.82, 'C', 'Color','b', 'FontSize', 22)
text(3.05, -0.15, 'B', 'Color','b', 'FontSize', 22)

% Add label for the fourth square
text(1.75, 2.15, 'A', 'Color','b', 'FontSize', 22)
text(2.9, 2.15, 'B', 'Color','b', 'FontSize', 22)

% Add label for the fifth square
text(3.9, 0.82, 'D', 'Color','b', 'FontSize', 22)
text(3.75, 2.15, 'A', 'Color','b', 'FontSize', 22)

% Add label for the sixth square
text(5.1, 1.85, 'B', 'Color','b', 'FontSize', 22)
text(4.9, 0.82, 'C', 'Color','b', 'FontSize', 22)

% Add label for the seventh square
text(4.9, 3.15, 'C', 'Color','b', 'FontSize', 22)
text(3.9, 3.15, 'D', 'Color','b', 'FontSize', 22)

% Add label for the eighth square
text(6, 3.15, 'D', 'Color','b', 'FontSize', 22)
text(6, 1.85, 'A', 'Color','b', 'FontSize', 22)

% Add label for the ninth square
%text(7.1, 3.15, 'C', 'Color','b', 'FontSize', 22)
%text(7.1, 1.85, 'B', 'Color','b', 'FontSize', 22)

% Add label for the tenth square
%text(7.1, 4.15, 'B', 'Color','b', 'FontSize', 22)
%text(5.85, 4.15, 'A', 'Color','b', 'FontSize', 22)

% Dashed line for triangle
plot([0.2 4.2], [0 0], 'k--', 'linewidth', 3) 
plot([4.2 4.2], [0 2], 'k--', 'linewidth', 3)

% Initial angle label
ang([0.2 0], 0.3, [0 atan(1/2)],'k-')
text(0.65, 0.15, '\alpha', 'Color', 'k', 'FontSize', 18)

% Position labels
text(0.1, 0.3, 'x_0', 'Color', 'k', 'FontSize', 16)
text(4.1, 2.3, 'x_0', 'Color', 'k', 'FontSize', 16)

% For dissertation
%text(1, 2, '$p$', 'Color','r', 'FontSize', 22,'Interpreter','latex')
text(1, 2, '$\alpha$', 'Color','k', 'FontSize', 22,'Interpreter','latex')

% Controlling the plot
set(gca, 'linewidth', 3, 'FontSize', 22);
set(gcf,'color','w');  % background color to white instead of usual grey
set(gca,'xtick', [], 'ytick', [], 'visible','off')  % removes axes and ticks



%% PRE-IMAGES images for dissertation

% Period-6 orbit with maximal area
init_pos = 0.25; init_angle = atan(2);
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'o-', false);
P_anim = square_billiards_draw(side, alpha, position); hold on

plot(init_pos, 0, 'ro', 'MarkerFaceColor', 'r')
text(-0.2, 1, '(a)', 'FontSize', 22)

%% Period-6 orbit that narrowly misses the vertices
init_pos = 0.4; init_angle = atan(2);
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'o-', false);
P_anim = square_billiards_draw(side, alpha, position); hold on

plot(init_pos, 0, 'ro', 'MarkerFaceColor', 'r')
text(-0.2, 1, '(b)', 'FontSize', 22)


%%
% Period-6 slope that hits vertices

init_pos = 2; init_angle = atan(2);
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'o-', false);
P_anim = square_billiards_draw(side, alpha, position); hold on

plot(0.5, 0, 'ro', 'MarkerFaceColor', 'r') 
text(-0.2, 1, '(c)', 'FontSize', 22)


%%
% Period-10 slope that hits vertices
init_pos = 2; init_angle = atan(2/3);
[side alpha position] = square_billiards_cobweb(init_angle, init_pos, N, 'o-', false);
P_anim = square_billiards_draw(side, alpha, position);

%% INFINTE TILING FIGURE FOR DISSERTATION. GRID. Figure 2.7

% Horizontal lines
plot([0 5], [0 0], 'b', 'linewidth', 3); hold on % 0<x<5 and 0<y<0
plot([0 5], [1 1], 'b', 'linewidth', 3)
plot([0 5], [2 2], 'b', 'linewidth', 3)
plot([0 5], [3 3], 'b', 'linewidth', 3)
plot([0 5], [4 4], 'b', 'linewidth', 3)
plot([0 5], [5 5], 'b', 'linewidth', 3)
% Vertical lines
plot([0 0], [0 5], 'b', 'linewidth', 3)
plot([1 1], [0 5], 'b', 'linewidth', 3)
plot([2 2], [0 5], 'b', 'linewidth', 3)
plot([3 3], [0 5], 'b', 'linewidth', 3)
plot([4 4], [0 5], 'b', 'linewidth', 3)
plot([5 5], [0 5], 'b', 'linewidth', 3)

% Sample trajectories
% Start from (0,0) and end up at (4,5)
x = [0:5]; y = 5/4 * x;
plot(x, y, 'r', 'linewidth', 2)

% Start from (0,0) and end up at (1,4)
x = [0:1]; y = 4 * x;
plot(x, y, 'g--', 'linewidth', 2)

% Start from (0,0) and end up at (3,2)
x = [0:0.1:3]; y = 2/3 * x;
plot(x, y, 'm.', 'linewidth', 2)

pbaspect([5 5 1])  % aspect ratio
set(gca,'XLim',[0 5], 'YLim', [0 5])

text(-0.25, -0.2, 'A', 'Color', 'b', 'FontSize', 22)
text(1.05, -0.2, 'B', 'Color','b', 'FontSize', 22)
text(1.05, 1.15, 'C', 'Color','b', 'FontSize', 22)
text(-0.25, 1.15, 'D', 'Color','b', 'FontSize', 22)

% For Proposition 2.4.4 where we cannot return with angle pi-\alpha
x = [0:0.1:4];
y = 1/(2-0.2)*(x-0.2);
plot(x, y, 'k', 'linewidth', 2)
plot(0.2, 0, 'ko', 'MarkerFaceColor', 'k')
plot(3.8, 2, 'ko', 'MarkerFaceColor', 'k')
text(0.2, -0.2, '$P_0$', 'FontSize', 18)
text(3.75, 1.8, '$P_0$', 'FontSize', 18)

text(4.05, 2.15, 'A', 'Color', 'b', 'FontSize', 22)cle
text(2.75, 2.15, 'B', 'Color','b', 'FontSize', 22)
text(2.75, 3.15, 'C', 'Color','b', 'FontSize', 22)
text(4.05, 3.15, 'D', 'Color','b', 'FontSize', 22)

ang([0.2 0], 0.4, [0 atan(0.5)],'k-');
text(0.65, 0.15, '$\alpha_0$', 'Color','k', 'FontSize', 18)

% Controlling the plot
set(gca, 'linewidth', 3, 'FontSize', 22, 'FontName', 'Helvetica');
set(gcf,'color','w');  % background color to white instead of usual grey
set(gca,'xtick', [], 'ytick', [], 'visible','off')  % removes axes and ticks



%% SHOULD WRITE FUNCTION IN ANOTHER FILE PURELY FOR THE ANIMATION

N = 60;
[side alpha position] = square_billiards_cobweb(2.6015, P_0, N, 'o', false);


P_anim = [position; zeros(size(position))]';   % matrix of x and y co-ords as columns

for i = 1:(N+1)
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

plot(P_anim(:,1), P_anim(:,2), 'k', 'linewidth', 1), hold on, %pbaspect([1 1 1])  % aspect ratio

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
set(gca,'XLim',[0 1], 'YLim', [0 1],'xtick', [], 'ytick', [])

% Formatting title
%th = title(sprintf('Trajectory for initial angle of %0.2f and initial position of %0.2f', init_angle, init_pos))
%titlePos = get( th , 'position');
%titlePos(2) = 1.05;
%set( th , 'position' , titlePos);

set(gcf,'color','w');  % background color to white instead of usual grey


%% SHOULD WRITE FUNCTION IN ANOTHER FILE PURELY FOR THE ANIMATION

N = 30;
[side alpha position] = square_billiards_cobweb(pi-atan(3/5), P_0, N, 'o', false);


P_anim = [position; zeros(size(position))]';   % matrix of x and y co-ords as columns

for i = 1:(N+1)
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

plot(P_anim(:,1), P_anim(:,2), 'k', 'linewidth', 1), hold on, %pbaspect([1 1 1])  % aspect ratio

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
set(gca,'XLim',[0 1], 'YLim', [0 1],'xtick', [], 'ytick', [])

% Formatting title
%th = title(sprintf('Trajectory for initial angle of %0.2f and initial position of %0.2f', init_angle, init_pos))
%titlePos = get( th , 'position');
%titlePos(2) = 1.05;
%set( th , 'position' , titlePos);

set(gcf,'color','w');  % background color to white instead of usual grey


