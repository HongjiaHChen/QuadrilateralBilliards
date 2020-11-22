%% RECTANGLE TESTING
clear all; clc

L_H = vpa(1); L_V = vpa(3.5);

x = linspace(0.01, L_H-0.01, 100);   % Initial positions on base side

period = 10;

for j=1:length(x)
    y(j) = rec_area(L_H, L_V, vpa(atan(4*L_V/1*L_H)), vpa(x(j)), period);
 
end

plot(x,y)
xlabel('P_0'); ylabel('Area')


%% Rectangle Newton solver

clear all; clc
eps = 1e-6;   % for finite differencing

% Size of rectangle
L_H = vpa(1); L_V = vpa(3.16);

period = 10; % what period are we seaching for?

alpha = 0.98; P = 0.49;

guesses = [alpha P]; % First column: init_angle and init_pos  

for j=1:7
    guess = guesses(j,:)';  % column vector
    
    % F^{N}(\alpha, \P) term. Where are we after a period number of bounces
    [F_side F_alpha, F_P] = rec_map(L_H, L_V, alpha, P, period);
    F_n = [F_alpha(period+1); F_P(period+1)];
    
    dA_dP = (rec_area(L_H, L_V, alpha, P+eps, period)-rec_area(L_H, L_V, alpha, P-eps, period))/(2*eps);

    guess = guess - inv(rec_jacobian_area(L_H, L_V, alpha, P, period)) * [F_n(2)-P; dA_dP];

    alpha = guess(1); P = guess(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TEMP SOLUTION for when we hop out of range
    alpha = mod(alpha, vpa(pi/2)); P = mod(P, vpa(1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    guesses(j+1,:) = [alpha P];
    
    if norm(guesses(j+1,:)- guesses(j,:))/norm(guesses(j,:)) < 1e-6   % relative error
        break
    end
end

guesses








%% FOLLOWING A PERIOD-4 ORBIT
clear all; clc
tic
period = 10; % what period are we seaching for?

% We know for square that the angle and position that is a period-4 orbit
% and that maximises the area is angle=pi/4, pos = 0.5. This is our
% initial guess.

eps = 1e-6;   % for finite differencing

height_perturb = 5e-2; % perturbing height of rectangle

heights = [1];             % The heights, filled in for square already
alpha_P_mat = [atan(1/4) 0.5]; % Already filled in for the square

areas = square_area_V2(alpha_P_mat(1), alpha_P_mat(2), period);   % store the areas as well
newton_steps = [0];                        % number of steps to convergence
function_eval = zeros(1,2);                       % what the function evaluates to (should be close to 0).

for i=1:50   % gradually change the height of the rectangle
    
    heights(i+1) = 1+height_perturb*i;
    L_H = vpa(1); L_V = vpa(heights(i+1)); % Size of rectangle

    guesses = alpha_P_mat(i,:); % First column: init_angle and init_pos  

    for j=1:10
        guess = guesses(j,:)';  % column vector
        alpha = guess(1); P = guess(2);

        % F^{N}(\alpha, \P) term. Where are we after a period number of bounces
        [F_side F_alpha, F_P] = rec_map(L_H, L_V, alpha, P, period);
        F_n = [F_alpha(period+1); F_P(period+1)];

        dA_dP = (rec_area(L_H, L_V, alpha, P+eps, period)-rec_area(L_H, L_V, alpha, P-eps, period))/(2*eps);

        guess = guess - inv(rec_jacobian_area(L_H, L_V, alpha, P, period)) * [F_n(2)-P; dA_dP];

        alpha = guess(1); P = guess(2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TEMP SOLUTION for when we hop out of range
        alpha = mod(alpha, vpa(pi/2)); P = mod(P, vpa(1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        guesses(j+1,:) = [alpha P];

        if (norm(guesses(j+1,:)- guesses(j,:))/norm(guesses(j,:)) < 1e-6)  ||  (j == 10)% relative error
            
            areas(i+1) = rec_area(L_H, L_V, alpha, P, period);
            newton_steps(i+1) = j;    % we converged on jth iteration
            
            % For function evaluation
            [F_side F_alpha, F_P] = rec_map(L_H, L_V, alpha, P, period);
            F_n = [F_alpha(period+1); F_P(period+1)];
            dA_dP = (rec_area(L_H, L_V, alpha, P+eps, period)-rec_area(L_H, L_V, alpha, P-eps, period))/(2*eps);
            function_eval(i+1,:) = [F_n(2)-guesses(j+1,2) dA_dP];
            
            break
        end
    end
    
    alpha_P_mat(i+1,:) = guesses(end,:);  % the converged value
end

bigmat = [heights' alpha_P_mat areas' newton_steps' function_eval];
toc


%% TEMP BACKWARDS
tic
period = 10; % what period are we seaching for?

% We know for square that the angle and position that is a period-4 orbit
% and that maximises the area is angle=pi/4, pos = 0.5. This is our
% initial guess.

eps = 1e-6;   % for finite differencing

height_perturb = 1e-2; % perturbing height of rectangle

heights_back = [1];             % The heights, filled in for square already
alpha_P_mat_back = [atan(1/4) 0.5]; % Already filled in for the square


areas_back = square_area_V2(alpha_P_mat_back(1), alpha_P_mat_back(2), period);   % store the areas as well
newton_steps_back = [0];                        % number of steps to convergence
function_eval_back = zeros(1,2);                       % what the function evaluates to (should be close to 0).


for i=1:30   % gradually change the height of the rectangle
    
    heights_back(i+1) = 1-height_perturb*i;
    L_H = vpa(1); L_V = vpa(heights_back(i+1)); % Size of rectangle

    guesses = alpha_P_mat_back(i,:); % First column: init_angle and init_pos  

    for j=1:10
        guess = guesses(j,:)';  % column vector
        alpha = guess(1); P = guess(2);

        % F^{N}(\alpha, \P) term. Where are we after a period number of bounces
        [F_side F_alpha, F_P] = rec_map(L_H, L_V, alpha, P, period);
        F_n = [F_alpha(period+1); F_P(period+1)];

        dA_dP = (rec_area(L_H, L_V, alpha, P+eps, period)-rec_area(L_H, L_V, alpha, P-eps, period))/(2*eps);

        guess = guess - inv(rec_jacobian_area(L_H, L_V, alpha, P, period)) * [F_n(2)-P; dA_dP];

        alpha = guess(1); P = guess(2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TEMP SOLUTION for when we hop out of range
        alpha = mod(alpha, vpa(pi/2)); P = mod(P, vpa(1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        guesses(j+1,:) = [alpha P];

        if norm(guesses(j+1,:)- guesses(j,:))/norm(guesses(j,:)) < 1e-6  ||  (j == 10)% relative error
            
            areas_back(i+1) = rec_area(L_H, L_V, alpha, P, period);
            newton_steps_back(i+1) = j;    % we converged on jth iteration
            
            % For function evaluation
            [F_side F_alpha, F_P] = rec_map(L_H, L_V, alpha, P, period);
            F_n = [F_alpha(period+1); F_P(period+1)];
            dA_dP = (rec_area(L_H, L_V, alpha, P+eps, period)-rec_area(L_H, L_V, alpha, P-eps, period))/(2*eps);
            function_eval_back(i+1,:) = [F_n(2)-guesses(j+1,2) dA_dP];
            
            break
        end
    end
    
    alpha_P_mat_back(i+1,:) = guesses(end,:);  % the converged value
end

bigmat_back = [heights_back' alpha_P_mat_back areas_back' newton_steps_back' function_eval_back];
toc


%% Combining the two matrices and saving

period10mat = [flip(bigmat_back); bigmat(2:end,:)];
save('Data\rec_period10_alpha1_4.mat', 'period10mat');












%% PLOTS

plot(combined_heights, combined_alpha_P_mat(:,1), 'o')

xlabel('Height'); ylabel('\alpha_0')
set(gca,'XLim',[min(combined_heights) max(combined_heights)], 'YLim', [0 pi/2])

yticks([0 pi/4 pi/2])
yticklabels({'0','\pi/4', '\pi/2'})
title(sprintf('Period %d', period))







%% Plot of height vs. P

plot(combined_heights, combined_alpha_P_mat(:,2), 'o')  % Second column is P
xlabel('Height'); ylabel('P_0')
set(gca,'XLim',[min(combined_heights) max(combined_heights)], 'YLim', [0 1])






%% Loading pre-computed results


p4_heights = matfile('Data\period4_heights.mat');
p4_heights = p4_heights.combined_heights;

p4_alpha_p_mat = matfile('Data\period4_alpha_P.mat');
p4_alpha_p_mat = p4_alpha_p_mat.combined_alpha_P_mat;

% Plotting
plot(p4_heights, p4_alpha_p_mat(:,1), 'o') % Angle 

xlabel('Height'); ylabel('\alpha_0')
set(gca,'XLim',[min(p4_heights) max(p4_heights)], 'YLim', [0 pi/2])
yticks([0 pi/4 pi/2])
yticklabels({'0','\pi/4', '\pi/2'})
title('Period 4')

%%
plot(p4_heights, p4_alpha_p_mat(:,2), 'o') % Position 


%% Period-6: atan(2)
p6_heights = matfile('Data\period6_heights.mat');
p6_heights = p6_heights.combined_heights;

p6_alpha_p_mat = matfile('Data\period6_alpha_P.mat');
p6_alpha_p_mat = p6_alpha_p_mat.combined_alpha_P_mat;

% Plotting
plot(p6_heights, p6_alpha_p_mat(:,1), 'o') % Angle 

xlabel('Height'); ylabel('\alpha_0')
set(gca,'XLim',[min(p6_heights) max(p6_heights)], 'YLim', [0 pi/2])
yticks([0 pi/4 pi/2])
yticklabels({'0','\pi/4', '\pi/2'})
title('Period 6')


%% Period-8: atan(3)

p8_mat = matfile('Data\rec_period8_alpha_3_1.mat');
p8_mat = p8_mat.period8mat;

% First column is height
% Column 2 is alpha
% Column 3 is P
% Column 4 is area
% Column 5 is Newton iterations
% Column 6&7 are function evaluations

plot3(p8_mat(:,1), p8_mat(:,2), p8_mat(:,3), 'o')
set(gca,'XLim',[0 max(p8_mat(:,1))], 'ZLim',[0 pi/2]);
xlabel('Height of rectangle')
ylabel('\alpha_0 - angle')
zlabel('P_0 - position')
title('Continuation Period-8')

%% Period-8 Area Plot

plot(p8_mat(:,1), p8_mat(:,4), 'o')
xlabel('Height of rectangle')
ylabel('Area of trajectory')
title('Area vs. height of rectangle')

%% Period-8 Newton iteration count plot

plot(p8_mat(:,1), p8_mat(:,5), 'o')
xlabel('Height of rectangle')
ylabel('Number of Newton iterations')
title('Number of Newton iterations vs. height of rectangle')

%% Period 8 Function evaluations

plot3(p8_mat(:,1), p8_mat(:,6), p8_mat(:,7), 'o')
set(gca,'XLim',[0 max(p8_mat(:,1))]);
xlabel('Height of rectangle')
ylabel('g^{N}(\alpha, P) - P')
zlabel('dA/dP')
title('Period-8 Function evaluations')





%% Period-10: For Dissertation


p10_alpha4_1_mat = matfile('Data\rec_period10_alpha4_1.mat');
p10_alpha3_2_mat = matfile('Data\rec_period10_alpha3_2.mat');
p10_alpha1_4_mat = matfile('Data\rec_period10_alpha1_4.mat');
p10_alpha2_3_mat = matfile('Data\rec_period10_alpha2_3.mat');

p10_alpha4_1_mat = p10_alpha4_1_mat.period10mat;
p10_alpha3_2_mat = p10_alpha3_2_mat.period10mat;
p10_alpha1_4_mat = p10_alpha1_4_mat.period10mat;
p10_alpha2_3_mat = p10_alpha2_3_mat.period10mat;

% alpha = atan(4), P = 0.125, 0.375, 0.625, 0.875
plot3(p10_alpha4_1_mat(:,1), p10_alpha4_1_mat(:,2), p10_alpha4_1_mat(:,3), 'b^', 'MarkerFaceColor', 'b'); hold on  % P_0 = 0.125
plot3(p10_alpha4_1_mat(:,1), atan(4.*p10_alpha4_1_mat(:,1)./1), repmat(0.375, size(p10_alpha4_1_mat(:,1))), 'y<', 'MarkerFaceColor', 'y'); hold on
plot3(p10_alpha4_1_mat(:,1), atan(4.*p10_alpha4_1_mat(:,1)./1), repmat(0.625, size(p10_alpha4_1_mat(:,1))), 'gv', 'MarkerFaceColor', 'g'); hold on
plot3(p10_alpha4_1_mat(:,1), atan(4.*p10_alpha4_1_mat(:,1)./1), repmat(0.875, size(p10_alpha4_1_mat(:,1))), 'r>', 'MarkerFaceColor', 'r'); hold on

% alpha = atan(3/2), P = 1/6, 0.5, 5/6
% plot3(p10_alpha3_2_mat(:,1), p10_alpha3_2_mat(:,2), p10_alpha3_2_mat(:,3), 'ro'); hold on
% plot3(p10_alpha3_2_mat(:,1), atan(3.*p10_alpha3_2_mat(:,1)./2), repmat(5/6, size(p10_alpha3_2_mat(:,1))), 'r*'); hold on
% plot3(p10_alpha3_2_mat(:,1), atan(3.*p10_alpha3_2_mat(:,1)./2), repmat(1/6, size(p10_alpha3_2_mat(:,1))), 'rs'); hold on

% alpha = atan(1/4), P = 0.5
%plot3(p10_alpha1_4_mat(:,1), p10_alpha1_4_mat(:,2), p10_alpha1_4_mat(:,3), 'gd'); hold on
x = linspace(0.01, 11, 100);
% plot3(x, atan(x./4), 0.5*ones(size(x)), 'gd'); hold on

indices = 1:1:size(p10_alpha2_3_mat,1);
% alpha = atan(2/3), P = 0.25,0.75
plot3(p10_alpha2_3_mat(indices,1), p10_alpha2_3_mat(indices,2), p10_alpha2_3_mat(indices,3), 'mo', 'MarkerFaceColor', 'm'); hold on
plot3(p10_alpha2_3_mat(:,1), atan(2.*p10_alpha2_3_mat(:,1)./3), repmat(3/4, size(p10_alpha2_3_mat(:,1))), 'ch', 'MarkerFaceColor', 'c'); hold on
%plot3(x, atan(2.*x./3), 0.75*ones(size(x)), 'mh'); hold on

set(gca,'XLim',[0 max(p10_alpha3_2_mat(:,1))], 'YLim',[0 pi/2], 'ZLim',[0 1], 'FontSize', 15);
xlabel('h - aspect ratio', 'FontSize', 15)
ylabel('\alpha_0 - angle', 'FontSize', 15)
zlabel('P_0 - position', 'FontSize', 15)
zticks([0, 0.25, 0.5, 0.75, 1])
zticklabels({'0', '0.25', '0.5', '0.75', '1'})
yticks([0 pi/4 pi/2])
yticklabels({'0','\pi/4', '\pi/2'})
title('Continuation Period-10')

grid on % background grid

% legend({'\alpha_0=atan(4A), P_0=0.125','\alpha_0=atan(4A), P_0=0.375', '\alpha_0=atan(4A), P_0=0.625','\alpha_0=atan(4A), P_0=0.875',...
%     '\alpha_0=atan(3A/2), P_0 = 1/6', '\alpha_0=atan(3A/2), P_0 = 0.5', '\alpha_0=atan(3A/2), P_0 = 5/6', ...
%     '\alpha_0=atan(A/4), P_0 = 0.5','\alpha_0=atan(2A/3),P_0 = 0.25','\alpha_0=atan(2A/3),P_0 = 0.75'})

[~, objh] = legend({'$(\alpha_0, P_0)=(\mathrm{atan}4h, 0.125)$','$(\alpha_0, P_0)=(\mathrm{atan}4h, 0.375)$', ...
    '$(\alpha_0, P_0)=(\mathrm{atan}4h, 0.625)$','$(\alpha_0, P_0)=(\mathrm{atan}4h, 0.875)$',...
     '$(\alpha_0, P_0)=(\mathrm{atan}\frac{2h}{3}, 0.25)$','$(\alpha_0, P_0)=(\mathrm{atan}\frac{2h}{3}, 0.75)$'},'Interpreter','latex')

view(-50, 13)
set(gcf,'color','w'); % change background color to white
legend boxoff   % hide box around the legend

objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
set(objhl, 'Markersize', 12); %// set marker size as desired


%%


%% Area                                                  alpha = atan(4)

plot(p10_alpha4_1_mat(:,1), p10_alpha4_1_mat(:,4), 'o')
xlabel('Aspect ratio of rectangle', 'FontSize', 15)
ylabel('Area of trajectory', 'FontSize', 15)
title('Area vs. aspect ratio of rectangle')


%% Newton steps                                          alpha = atan(4)

plot(p10_alpha4_1_mat(:,1), p10_alpha4_1_mat(:,5), 'o')
xlabel('Aspect ratio of rectangle', 'FontSize', 15)
ylabel('Number of Newton steps', 'FontSize', 15)
title('Number of Newton steps vs. aspect ratio of rectangle')
yticks([0 1 2 3 4 5])
yticklabels({'0','1','2','3','4','5'})
set(gca, 'YLim',[0 max(p10_alpha4_1_mat(:,5))+1]);

%% Function evaluation                                        alpha = atan(4)

plot3(p10_alpha4_1_mat(:,1), p10_alpha4_1_mat(:,6), p10_alpha4_1_mat(:,7), 'o')
xlabel('Aspect ratio of rectangle', 'FontSize', 15)
ylabel('$g^{(N)}(\alpha,P)-P$', 'FontSize', 15,'Interpreter','latex')
zlabel('$\frac{dA}{dP}$', 'FontSize', 18,'Interpreter','latex')
title('Function evaluation vs. aspect ratio of rectangle')







%%

plot3(p10_heights, p10_alpha_p_mat(:,1), p10_alpha_p_mat(:,2), 'o')
set(gca,'XLim',[0 max(p10_heights)], 'YLim',[0 pi/2]);

xlabel('Height of rectangle')
%xticks([0 1 2 3 4])
%xticklabels({'0', '1', '2', '3', '4'})

ylabel('\alpha_0 - angle')

zlabel('P_0 - position')







%% Angle of atan(1/4)

p10_alpha1_4_mat = matfile('Data\rec_period10_alpha1_4.mat');
p10_alpha1_4_mat = p10_alpha1_4_mat.period10mat;

plot3(p10_alpha1_4_mat(:,1), p10_alpha1_4_mat(:,2), p10_alpha1_4_mat(:,3), 'o'); hold on
set(gca,'XLim',[0 max(p10_alpha1_4_mat(:,1))], 'YLim',[0 pi/2], 'ZLim',[0 1]);
xlabel('Height of rectangle')
ylabel('\alpha_0 - angle')
zlabel('P_0 - position')
title('Continuation Period-10')


plot3(p10_alpha1_4_mat(:,1), atan(1.*p10_alpha1_4_mat(:,1)./4), p10_alpha1_4_mat(:,3), 'r-')


%% Area
plot(p10_alpha1_4_mat(:,1), p10_alpha1_4_mat(:,4), 'o')
xlabel('Height of rectangle')
ylabel('Area of trajectory')
title('Area vs. height of rectangle')

%% newton steps
plot(p10_alpha1_4_mat(:,1), p10_alpha1_4_mat(:,5), 'o')
xlabel('Height of rectangle')
ylabel('Area of trajectory')
title('Area vs. height of rectangle')
