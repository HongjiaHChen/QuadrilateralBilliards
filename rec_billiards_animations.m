clear all; close all; clc; 

% Animations of the billiard trajectory

% Our Difference Equation for Rectangular Billiards:
% x(n+1) = f(x(n)). f is too complex to write here. Refer to written notes.

N = 100;   % number of iterations/maps (bounces)

% NEED TO USE vpa() FOR SYMBOLIC ACCURACY

%% TESTING DELETE






L_H = vpa(1); L_V = vpa(0.5); % Horizontal (base) length and Vertical (height) length
init_pos = vpa(0.5); init_angle = vpa(atan(L_V/L_H));%init_angle = vpa(pi/4); 
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);









%% Rectangle with base length of 2 and height of 1.

L_H = vpa(2); L_V = vpa(1); % Horizontal (base) length and Vertical (height) length
init_pos = vpa(0.5); init_angle = vpa(atan(L_V/L_H));%init_angle = vpa(pi/4); 
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);

%% Rectangle with base length of 4 and height of 3.
% Period 10 orbit
L_H = vpa(4); L_V = vpa(3); % Horizontal (base) length and Vertical (height) length
init_pos = vpa(0.5); init_angle = vpa(atan(1/2));%init_angle = vpa(pi/4); 
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);


%% Rectangle with base length of 1 and height of 2.2.
% Period 32 orbit with lines and dots. Initial angle of pi/4 and position of 0.7. 

L_H = vpa(1); L_V = vpa(2.2); % Horizontal (base) length and Vertical (height) length
init_pos = vpa(0.7); init_angle = vpa(pi/4);
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);


%% Rectangle with base length of 2.2 and height of 3.5.
% Period 6 orbit
L_H = vpa(2.2); L_V = vpa(3.5); % Horizontal (base) length and Vertical (height) length
init_pos = vpa(0.5); init_angle = atan(vpa((2*3.5)/(3*2.2)));  % numerator: 2*base. denominator: 2*height
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);


%% Keeping same initial position 0.5 and initial angle pi/4. And perturbing the height.

init_pos = sym(0.5); init_angle = sym(pi/4); L_H = sym(1);% KEEP THIS FIXED

L_V = sym(1); % THIS WILL CHANGE
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);

%% Height of 1.01
L_V = sym(1.01); % THIS WILL CHANGE
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);

%% Height of 1.01234
L_V = sym(1.01234); % THIS WILL CHANGE
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);

%% Height of 1.012345. THIS GIVES THE MOD ERRORS WITH THE ADDITIONAL DP. NOT SURE WHY SYMBOLIC VPA ADDS ON EXTRA DIGITS
% for L_V.
% Eg: check vpa(1.01234) compared to vpa(1.012345)
L_V = sym(1.012345); % THIS WILL CHANGE
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);

%% Height of 1.0324
L_V = vpa(1.0324); % THIS WILL CHANGE
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);

%% Height of 1.05
L_V = vpa(1.05); % THIS WILL CHANGE
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);

%% Height of 1.0745214
L_V = 1.0745214; % THIS WILL CHANGE
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);

%% Height of 1.1
L_V = 1.1; % THIS WILL CHANGE
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);

%% Height of 1.15
L_V = 1.15; % THIS WILL CHANGE
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);

%% Height of 1.499999.      THIS IS APPARENTLY PERIOD 10
L_V = 1.499999; % THIS WILL CHANGE
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);

%% Height of 1.5
L_V = 1.5; % THIS WILL CHANGE
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);

%% Height of 1.50001
L_V = 1.50001; % THIS WILL CHANGE
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);

%% Height of 1.51
L_V = 1.51; % THIS WILL CHANGE
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);


%%
L_V = sym(heights_vec(667)); % THIS WILL CHANGE
[side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position);







%% Loop through and produce images like Frank's

init_pos = vpa(0.5); init_angle = vpa(pi/4); L_H = vpa(1);% FIXED
N = 50; % how long the billiard goes for

N_iter = 2000;
heights_vec = linspace(0.5, 1.5, N_iter);

for j=1:N_iter  % change the heights with each iteration
    L_V = heights_vec(j); 
    [side alpha position] = rec_billiards_cobweb(L_H, L_V, init_angle, init_pos, N, 'o-', false);
    P_anim = rec_billiards_draw(L_H, L_V, side, alpha, position); hold off
    
    Image = getframe(gcf);
    folder = 'C:/Users/Work/Desktop/UoA/Year 4/Honours Project/Mathematical Billiards/MATLAB Code/Images/Square_Rectangle Bifurcations/Animation/50 iterations init_pos=0.5 init_angle = pi.4';
    baseFileName = sprintf('test%d.jpg', j);
    fullFileName = fullfile(folder, baseFileName);
    imwrite(Image.cdata, fullFileName);
end




%% Rectangle tesselation (for dissertation)

% Set a = 2, b = 1. Have a as the base of the rectangle and b has the
% height.
a = 2; b = 1;

% Horizontal lines
plot([0 5*a], [0 0], 'b', 'linewidth', 3); hold on  % 0<x<10, y=0 fixed
plot([0 5*a], [b b], 'b', 'linewidth', 3)           % Adding 0.5 to y co-ord
plot([0 5*a], [2*b 2*b], 'b', 'linewidth', 3) 
plot([0 5*a], [3*b 3*b], 'b', 'linewidth', 3) 
%plot([0 5*a], [4*b 4*b], 'b', 'linewidth', 3)  
%plot([0 5*a], [5*b 5*b], 'b', 'linewidth', 3)

% Vertical lines
plot([0 0], [0 5*b], 'b', 'linewidth', 3)           % x=0 fixed, 0<y<2.5
plot([a a], [0 5*b], 'b', 'linewidth', 3)
plot([2*a 2*a], [0 5*b], 'b', 'linewidth', 3)
plot([3*a 3*a], [0 5*b], 'b', 'linewidth', 3)
plot([4*a 4*a], [0 5*b], 'b', 'linewidth', 3)
plot([5*a 5*a], [0 5*b], 'b', 'linewidth', 3)

% text(1, 2.5, '$b$', 'Color','g', 'FontSize', 22,'Interpreter','latex') % for diss

% Straightline trajectory
x = [0:5*a]; y = 0.25*x - 0.1;
plot(x, y, 'k', 'linewidth', 2)

% Dashed line for triangle
plot([0.4 8.4], [0 0], 'k--', 'linewidth', 3) 
plot([8.4 8.4], [0 2], 'k--', 'linewidth', 3)

% Initial angle label
ang([0.4 0], 0.9, [0 atan(1/4)],'k-')
text(1.5, 0.2, '\alpha', 'Color', 'k', 'FontSize', 18)

% Position labels
text(0.2, 0.3, 'x_0', 'Color', 'k', 'FontSize', 16)
text(8.1, 2.3, 'x_0', 'Color', 'k', 'FontSize', 16)

% Size of rectangle labels
text(1, 1.4, '$a$', 'Color', 'g', 'FontSize', 22,'Interpreter','latex')
text(-0.6, 0.5, '$b$', 'Color', 'g', 'FontSize', 22,'Interpreter','latex')

% Number of flips labels
text(5.5, -0.8, '$q$ flips', 'Color','r', 'FontSize', 22,'Interpreter','latex')
h = text(-0.6, 1.6, '$p$ flips', 'Color','r', 'FontSize', 22,'Interpreter','latex')
set(h,'Rotation',90);

pbaspect([5*a 3*b 1])  % aspect ratio
set(gca,'XLim',[0 5*a], 'YLim', [0 3*b])

% Controlling the plot
set(gca, 'linewidth', 3, 'FontSize', 22);
set(gcf,'color','w');  % background color to white instead of usual grey
set(gca,'xtick', [], 'ytick', [], 'visible','off')  % removes axes and ticks

