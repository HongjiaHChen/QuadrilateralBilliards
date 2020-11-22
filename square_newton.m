%% Given a particular period and initial angle. We want the initial position which maximises
% the 'area' of the trajectory. The definition of the area is ambiguous and
% care must be taken in choosing it.

% We use the definition in square_area_V2.m. We find that the pre-images of
% the vertex and vertices give minimal 'area' with this definiton and
% midpoints of pre-images and vertices give maximal area.

period = 10; initial_angle = atan(4/1);

x = linspace(0.001, 0.999, 10000);

for j=1:length(x)
    y(j) = square_area_V2(initial_angle, x(j), period);
end

plot(x,y,'LineWidth',4)
xlabel('Initial Position - P_0', 'FontSize', 15); ylabel('Area', 'FontSize', 15)
set(gca,'FontSize',20)



%% Derivative of our Area function using Central difference method

alpha_star = pi/4; P_star = 0.4; period = 4;

eps = 1e-6;   % for finite differencing

A_big_P = square_area_V2(alpha_star, P_star+eps, period);
A_small_P = square_area_V2(alpha_star, P_star-eps, period);

dA_dP = (A_big_P-A_small_P)/(2*eps)



%% Finding alpha and P which maximise the area and have that P = g^(N)(alpha,P). i.e: points which 
% return to P after N iterations.

clear all; clc
eps = 1e-6;   % for finite differencing

period = 10; % what period are we seaching for?

alpha = 1.25; P = 0.24;
guesses = [alpha P]; % First column: init_angle and init_pos   


for j=1:7   % Only allow a maximum of 7 iterations
    guess = guesses(j,:)';  % column vector
    
    % F^{N}(\alpha, \P) term.Where are we after a period number of bounces?
    [F_alpha, F_P] = square_map(alpha, P, period);
    F_n = [F_alpha(period+1); F_P(period+1)];
    
    dA_dP = (square_area_V2(alpha, P+eps, period)-square_area_V2(alpha, P-eps, period))/(2*eps);

    guess = guess - inv(square_jacobian_area(alpha, P, period)) * [F_n(2)-P; dA_dP]
    
    % det(square_jacobian_area(alpha, P, period))
    
    alpha = guess(1); P = guess(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TEMP SOLUTION for when we hop out of range
    alpha = mod(alpha, pi/2); P = mod(P, 1);
    guesses(j+1,:) = [alpha P];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Checking if our solution is accurate enough
    if norm(guesses(j+1,:)- guesses(j,:))/norm(guesses(j,:)) < 1e-6   % relative error
        break
    end
end

guesses





%% Alternative approach. (NOT USED ANYMORE). Only treating P as a variable.
% We have the scalar equation: F(P) = P-g^(N)(alpha,P) + dA/dP(alpha,P)=0

% The P which solves F(P) maximises the area with the given alpha and N.

N = 8; % period
P = 0.1; alpha = atan(3);

h = 1e-6;

for j=1:3
    [N_alpha, N_P] = square_map(alpha, P+h, N); g_big = N_P(N+1);
    [N_alpha, N_P] = square_map(alpha, P-h, N); g_small = N_P(N+1);
    [N_alpha, N_P] = square_map(alpha, P, N);   g = N_P(N+1);

    A_big = square_area(alpha, P+h, N);
    A = square_area(alpha, P, N);
    A_small = square_area(alpha, P-h, N);

    F_prime  = 1 - (g_big-g_small)/(2*h) + (A_big-2*A+A_small)/(h^2);

    F_p = P-g+(A_big-A_small)/(2*h);
    
    
    P = P - F_p/F_prime
end




%% For the dissertation
% Period 4 area. alpha = atan(1)

period = 4; initial_angle = atan(1/1);

x = linspace(0.001, 0.999, 100000);

for j=1:length(x)
    y(j) = square_area_V2(initial_angle, x(j), period);
end

plot(x,y,'LineWidth',4)
xlabel('Initial Position - P_0', 'FontSize', 15); ylabel('Area', 'FontSize', 15)
set(gca,'FontSize',20) % change size of axis numbers
set(gcf,'color','w'); % change background color to white
xticks([0, 0.25, 0.5, 0.75, 1])
xticklabels({'0', '0.25', '0.5', '0.75', '1'})
yticks([3, 3.25, 3.5])
yticklabels({'3', '3.25', '3.5'})

text(min(x)-0.13, max(y), '(a)','HorizontalAlignment', 'right', 'Fontsize', 20)
pbaspect([1 1 1])

% Image = getimage(gcf);
% folder = 'C:/Users/Work/Desktop/UoA/Year 4/Honours Project/Mathematical Billiards/Dissertation/Images/Chapter 4';
% baseFileName = sprintf('area_period4.jpg', j);
% fullFileName = fullfile(folder, baseFileName);
% imwrite(Image.cdata, fullFileName);

%% Period 10 area. alpha = atan(4)

clear all; clc
period = 10; initial_angle = atan(4/1);

x = linspace(0.001, 0.999, 10000);

for j=1:length(x)
    y(j) = square_area_V2(initial_angle, x(j), period);
end

plot(x,y,'LineWidth',4)
xlabel('Initial Position - P_0', 'FontSize', 15); ylabel('Area', 'FontSize', 15)
set(gca,'FontSize',20) % change size of axis numbers
set(gcf,'color','w'); % change background color to white
xticks([0, 0.25, 0.5, 0.75, 1])
xticklabels({'0', '0.25', '0.5', '0.75', '1'})
yticks([9.74, 9.81, 9.88])
yticklabels({'9.74', '9.81', '9.88'})

text(min(x)-0.1, 9.88, '(b)','HorizontalAlignment', 'right', 'Fontsize', 20)

% Image = getimage(gcf);
% folder = 'C:/Users/Work/Desktop/UoA/Year 4/Honours Project/Mathematical Billiards/Dissertation/Images/Chapter 4';
% baseFileName = sprintf('area_period10.jpg', j);
% fullFileName = fullfile(folder, baseFileName);
% imwrite(Image.cdata, fullFileName);


