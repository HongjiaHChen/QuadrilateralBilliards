%% Area testing

h = 1; gamma = pi/2;

x = linspace(0.01, 1-0.01, 1000);   % Initial positions on base side

period = 8;
alpha = atan(1/3);     %atan(h/(1+h/tan(gamma)));

for j=1:length(x)
    y(j) = parallelogram_area(h, gamma, alpha, x(j), period);
 
end

plot(x,y)
xlabel('P_0'); ylabel('Area')


%% Maximum area testing
clear all; clc

h = 3.5; gamma = pi/2-0.01;

alpha0 = atan(h*3/2);
period = 10;


x = linspace(0.01, 1-0.01, 100);   % Initial positions on base side

for j=1:length(x)
    y(j) = parallelogram_area(h, gamma, alpha0, x(j), period);
end

plot(x,y)
xlabel('P_0'); ylabel('Area')

%%

%% Parallelogram Newton solver

clear all; clc
eps = 1e-6;   % for finite differencing

% Size of parallelogram
h = 3.16; gamma = pi/2;

period = 10; % what period are we seaching for?

alpha = 0.98; P = 0.49;   % Initial guess

guesses = [alpha P]; % First column: init_angle and init_pos  

for j=1:7
    guess = guesses(j,:)';  % column vector
    
    % F^{N}(\alpha, \P) term. Where are we after a period number of bounces
    [F_side F_alpha, F_P] = parallelogram_map(h, gamma, alpha, P, period);
    F_n = [F_alpha(period+1); F_P(period+1)];
    
    dA_dP = (parallelogram_area(h, gamma, alpha, P+eps, period)-parallelogram_area(h, gamma, alpha, P-eps, period))/(2*eps);

    guess = guess - inv(parallelogram_jacobian_area(h, gamma, alpha, P, period)) * [F_n(2)-P; dA_dP];

    alpha = guess(1); P = guess(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TEMP SOLUTION for when we hop out of range
    alpha = mod(alpha, pi); P = mod(P, 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    guesses(j+1,:) = [alpha P];
    
    if norm(guesses(j+1,:)- guesses(j,:))/norm(guesses(j,:)) < 1e-6   % relative error
        break
    end
end

guesses




%% Following a period-10 orbit by keeping gamma constant and varying the height.

clear all; clc
tic
period = 10; % what period are we seaching for?

gamma = pi/2; % keep this fixed for now

% We know for square that the angle and position that is a period-4 orbit
% and that maximises the area is angle=pi/4, pos = 0.5. This is our
% initial guess.

eps = 1e-6;   % for finite differencing

height_perturb = 1e-2; % perturbing height of rectangle

heights = [1];             % The heights, filled in for square already
alpha_P_mat = [atan(2/3) 0.125]; % Already filled in for the square                      % CHANGE HERE

areas = square_area_V2(alpha_P_mat(1), alpha_P_mat(2), period);   % store the areas as well
newton_steps = [0];                        % number of steps to convergence
function_eval = zeros(1,2);                       % what the function evaluates to (should be close to 0).

for i=1:1000   % gradually change the height of the rectangle
    
    heights(i+1) = 1+height_perturb*i;
    
    h = heights(i+1); % Height of parallelogram

    guesses = alpha_P_mat(i,:); % First column: init_angle and init_pos  

    for j=1:10
        guess = guesses(j,:)';  % column vector
        alpha = guess(1); P = guess(2);

        % F^{N}(\alpha, \P) term. Where are we after a period number of bounces
        [F_side F_alpha, F_P] = parallelogram_map(h, gamma, alpha, P, period);
        F_n = [F_alpha(period+1); F_P(period+1)];

        dA_dP = (parallelogram_area(h, gamma, alpha, P+eps, period)-parallelogram_area(h, gamma, alpha, P-eps, period))/(2*eps);

        guess = guess - inv(parallelogram_jacobian_area(h, gamma, alpha, P, period)) * [F_n(2)-P; dA_dP];

        alpha = guess(1); P = guess(2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TEMP SOLUTION for when we hop out of range
        alpha = mod(alpha, pi); P = mod(P, 1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        guesses(j+1,:) = [alpha P];

        if (norm(guesses(j+1,:)- guesses(j,:))/norm(guesses(j,:)) < 1e-6)  ||  (j == 10)% relative error
            
            % For function evaluation. Are we at [0,0] yet?
            [F_side F_alpha, F_P] = parallelogram_map(h, gamma, alpha, P, period);
            F_n = [F_alpha(period+1); F_P(period+1)];
            dA_dP = (parallelogram_area(h, gamma, alpha, P+eps, period)-parallelogram_area(h, gamma, alpha, P-eps, period))/(2*eps);
            
            if norm([F_n(2)-guesses(j+1,2) dA_dP]) < 1e-6
                areas(i+1) = parallelogram_area(h, gamma, alpha, P, period);
                newton_steps(i+1) = j;    % we converged on jth iteration
                function_eval(i+1,:) = [F_n(2)-guesses(j+1,2) dA_dP];
                break
            end
        end
    end
    
    alpha_P_mat(i+1,:) = guesses(end,:);  % the converged value
end

bigmat = [heights' alpha_P_mat areas' newton_steps' function_eval];
toc


%% Going BACKWARDS in height
tic

gamma = pi/2; % angle for parallelogram

period = 10; % what period are we seaching for?

eps = 1e-6;   % for finite differencing

height_perturb = 1e-3; % perturbing height of rectangle

heights_back = [1];             % The heights, filled in for square already
alpha_P_mat_back = [atan(2/3) 0.125]; % Already filled in for the square            % CHANGE HERE


areas_back = square_area_V2(alpha_P_mat_back(1), alpha_P_mat_back(2), period);   % store the areas as well
newton_steps_back = [0];                        % number of steps to convergence
function_eval_back = zeros(1,2);                       % what the function evaluates to (should be close to 0).


for i=1:990   % gradually change the height of the rectangle
    
    heights_back(i+1) = 1-height_perturb*i;
    
    h = heights_back(i+1);  % height of parallelogram

    guesses = alpha_P_mat_back(i,:); % First column: init_angle and init_pos  

    for j=1:10
        guess = guesses(j,:)';  % column vector
        alpha = guess(1); P = guess(2);

        % F^{N}(\alpha, \P) term. Where are we after a period number of bounces
        [F_side F_alpha, F_P] = parallelogram_map(h, gamma, alpha, P, period);
        F_n = [F_alpha(period+1); F_P(period+1)];

        dA_dP = (parallelogram_area(h, gamma, alpha, P+eps, period)-parallelogram_area(h, gamma, alpha, P-eps, period))/(2*eps);

        guess = guess - inv(parallelogram_jacobian_area(h, gamma, alpha, P, period)) * [F_n(2)-P; dA_dP];

        alpha = guess(1); P = guess(2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TEMP SOLUTION for when we hop out of range
        alpha = mod(alpha, pi/2); P = mod(P, 1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        guesses(j+1,:) = [alpha P];

        if norm(guesses(j+1,:)- guesses(j,:))/norm(guesses(j,:)) < 1e-6  ||  (j == 10)% relative error
            
            % For function evaluation
            [F_side F_alpha, F_P] = parallelogram_map(h, gamma, alpha, P, period);
            F_n = [F_alpha(period+1); F_P(period+1)];
            dA_dP = (parallelogram_area(h, gamma, alpha, P+eps, period)-parallelogram_area(h, gamma, alpha, P-eps, period))/(2*eps);
            
            if norm([F_n(2)-guesses(j+1,2) dA_dP]) < 1e-6
                areas_back(i+1) = parallelogram_area(h, gamma, alpha, P, period);
                newton_steps_back(i+1) = j;    % we converged on jth iteration
                function_eval_back(i+1,:) = [F_n(2)-guesses(j+1,2) dA_dP];
                break
            end
        end
    end
    
    alpha_P_mat_back(i+1,:) = guesses(end,:);  % the converged value
end

bigmat_back = [heights_back' alpha_P_mat_back areas_back' newton_steps_back' function_eval_back];
toc

%%

period10mat = [flip(bigmat_back); bigmat(2:end,:)];
%save('Data\rec_period10_alpha2_3.mat', 'period10mat');





%% This matches the result from rec_newton
plot3(period10mat(:,1), period10mat(:,2), period10mat(:,3), 'o'); hold on
set(gca,'XLim',[0 max(period10mat(:,1))], 'YLim',[0 pi/2], 'ZLim',[0 1]);
xlabel('Height of rectangle')
ylabel('\alpha_0 - angle')
zlabel('P_0 - position')
yticks([0 pi/4 pi/2])
yticklabels({'0','\pi/4', '\pi/2'})
title('Continuation Period-10')

plot3(period10mat(:,1), atan(4.*period10mat(:,1)./1), period10mat(:,3), 'r-')

%% Area

plot(period10mat(:,1), period10mat(:,4), 'o')
xlabel('Height of Parallelogram')
ylabel('Area of trajectory')
title('Area vs. height of rectangle')

%% Newton steps

plot(period10mat(:,1), period10mat(:,5), 'o')
xlabel('Height of Parallelogram')
ylabel('Number of Newton steps')
title('Number of Newton steps vs. height of Parallelogram')
yticks([0 1 2 3 4 5])
yticklabels({'0','1','2','3','4','5'})
set(gca, 'YLim',[0 max(period10mat(:,5))+1]);




%% Following a period-6 orbit by keeping h constant and varying gamma.
clear all; clc
tic
period = 6*1; % what period are we seaching for?

h = 1; % keep this fixed for now

% We know for square that the angle and position that is a period-6 orbit
% and that maximises the area is angle=atan(2), pos = 0.25. This is our
% initial guess.

gamma = pi/2;             % for the square
alpha = atan(2); P = 0.25;

test_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);
toc

% Column 1 is gamma
% Column 2 is initial angle, alpha
% Column 3 is initial position, P
% Column 4 is the area
% Column 5 is the number of Newton steps
% Column 6 & 7 are function evaluations for g^(N)(alpha,P)-P and dA/dP
% Column 8 is the function evaluation for alpha, f^(N)(alpha,P)-alpha
test_mat = test_mat(1:(end-1),:);

%%

index = size(test_mat,1);   % last row

gamma = test_mat(index,1);    % the last gamma that we converged, this is where we start from now
alpha = test_mat(index,2); P = test_mat(index,3);  % the angle and position in for which we last converged

period = 36;

temp_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);

temp_mat
        
%%
N = 50; % number of animation iterations

row = size(temp_mat,1) - 2;
h = 1; gamma = temp_mat(row, 1) ;

alpha_star = temp_mat(row, 2); P_star = temp_mat(row, 3);   
[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star+0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);

%%

test_mat = [test_mat; temp_mat(2:(end-1),:)];


%% The last continuation before non-convergence

N = 24; % number of animation iterations

row = size(test_mat,1) - 3;
h = 1; gamma = test_mat(row, 1) ;

alpha_star = test_mat(row, 2); P_star = test_mat(row, 3);   
[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star+0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);



%% Animating the continuation of orbits

for j=1:100:size(test_mat,1)  % iterating over each row in the matrix
    if ~isnan(test_mat(j, 2))
        N = 50; % iterations for animations

        h = 1; gamma = test_mat(j, 1);

        alpha_star = test_mat(j, 2); P_star = test_mat(j, 3);   

        [side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
        P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

        Image = getframe(gcf);
%         folder = 'C:/Users/Work/Desktop/UoA/Year 4/Honours Project/Mathematical Billiards/MATLAB Code/Images/Parallelogram/Animation/Period 6';
%         baseFileName = sprintf('test%d.jpg', j);
%         fullFileName = fullfile(folder, baseFileName);
%         imwrite(Image.cdata, fullFileName);
    end
end


%% To continue the periodic orbit once it passes through the vertex. We increase the period we want to follow.
% As a period-2 orbit is also a period-6 orbit. We find that at each
% bifurcaton for this periodic orbit, the period increases by 4.

h = 1; gamma = pi/2;    % dimensions of initial parallelogram, we keep h fixed

alpha = atan(2); P = 0.25;  % the angle and position in the square for which we start following from
% (note this normally gives least period-2 orbits in the square)
    
Nmax = 36;

for i=6:6:Nmax
    
    period = i;          % following a period-i orbit
    
    if i == 6   % first iteration
        test_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);
        test_mat = test_mat(1:(end-1),:);   % delete last row
        
        % to make all the areas comparable as we have our area func depends on period
        cum_area = test_mat(:,4) - period .* h; % the cumulative area. Area = period*h - cum
        test_mat(:,4) = prod(6:4:Nmax).*h - cum_area .* prod(6:4:Nmax)./period;  
        
        test_mat(:,4) = test_mat(:,4)./(prod(6:4:Nmax).*h)   % normlization, to make numbers smaller
        
        
        
    else % use the previous iterations results to start from
        index = size(test_mat,1);   % last row

        gamma = test_mat(index,1);    % the last gamma that we converged, this is where we start from now
        alpha = test_mat(index,2); P = test_mat(index,3);  % the angle and position in for which we last converged

        temp_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);
        
        %test_mat(:,9) = test_mat(:,4);
        
        % to make all the areas comparable as we have our area func depends on period
        cum_area = temp_mat(:,4) - period .* h; % the cumulative area. Area = period*h - cum
        temp_mat(:,4) = prod(6:4:Nmax).*h - cum_area .* prod(6:4:Nmax)./period;  
        temp_mat(:,4) = temp_mat(:,4)./(prod(6:4:Nmax).*h);   % normlization, to make numbers smaller
        
        test_mat = [test_mat; temp_mat(2:(end-1),:)];  % appending the new results to the matrix

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TEMP PLOTTING

    N = 100; % iterations for animations

    index = size(test_mat,1);
    gamma = test_mat(index, 1);
    alpha_star = test_mat(index, 2);
    P_star = test_mat(index, 3);

    [side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star-0.00, N);
    P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off
    w = waitforbuttonpress;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end


test_mat










%% 3D Plot of gamma, alpha_0, position0
plot3(test_mat(:,1), test_mat(:,2), test_mat(:,3), 'go'); hold on
set(gca,'XLim',[0 pi/2], 'YLim',[0 pi], 'ZLim',[0 2+2*h/max(test_mat(:,1))]);
xlabel('Gamma of parallelogram')
ylabel('\alpha_0 - angle')
zlabel('P_0 - position')
title(sprintf('Continuation Period-%d', period))


%%

%line curve of vertice position
gammas = linspace(0.01,pi/2, 500);
plot3(gammas, 1.3.*ones(size(gammas)), 1 + 1./sin(gammas), 'k-'); hold on


%% Area

plot(test_mat(:,1), test_mat(:,4), 'o')
xlabel('gamma of parallelogram')
ylabel('Area of trajectory')
title('Area vs. gamma of parallelogram')


%% Newton steps

plot(test_mat(:,1), test_mat(:,5), 'o')
xlabel('gamma of parallelogram')
ylabel('Number of Newton steps')
title('Number of Newton steps vs. gamma of parallelogram')








%% ANOTHER PERIOD-6

clear all; clc
period = 6*4; % what period are we seaching for?

h = 1; % keep this fixed for now

gamma = pi/2;             % for the square
alpha = atan(2); P = 0.75;

test_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);
test_mat = test_mat(1:(end-1),:)



%%

index = size(test_mat,1);   % last row

gamma = test_mat(index,1);    % the last gamma that we converged, this is where we start from now
alpha = test_mat(index,2); P = test_mat(index,3);  % the angle and position in for which we last converged

period = 14;

temp_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);

temp_mat

%test_mat = [test_mat; temp_mat(2:(end-1),:)];
        
%%
N = 50; % number of animation iterations

row = size(temp_mat,1) - 20;
h = 1; gamma = temp_mat(row, 1) ;

alpha_star = temp_mat(row, 2); P_star = temp_mat(row, 3);   
[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star+0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);

%%

test_mat = [test_mat; temp_mat(2:(end-1),:)];


%% The last continuation before non-convergence

N = 24; % number of animation iterations

row = size(test_mat,1) - 3;
h = 1; gamma = test_mat(row, 1) ;

alpha_star = test_mat(row, 2); P_star = test_mat(row, 3);   
[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star+0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);



%% Animating the continuation of orbits

for j=1:20:size(test_mat,1)  % iterating over each row in the matrix
    if ~isnan(test_mat(j, 2))
        N = 50; % iterations for animations

        h = 1; gamma = test_mat(j, 1);

        alpha_star = test_mat(j, 2); P_star = test_mat(j, 3);   

        [side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
        P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

        Image = getframe(gcf);
%         folder = 'C:/Users/Work/Desktop/UoA/Year 4/Honours Project/Mathematical Billiards/MATLAB Code/Images/Parallelogram/Animation/Period 6';
%         baseFileName = sprintf('test%d.jpg', j);
%         fullFileName = fullfile(folder, baseFileName);
%         imwrite(Image.cdata, fullFileName);
    end
end











%% Period-2 orbit bifurcating into period-6 orbit.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
period = 2* 6;  %2*6*10;  % following a period-6 orbit
h = 1; gamma = pi/2;    % dimensions of parallelogram

alpha = pi/2; P = 3.5;  % the angle and position in the square for which we start following from
% (note this normally gives least period-2 orbits in the square)

test_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);

test_mat

%% The animation for the period-2 to period-6 bifurcation

for j=1:10:size(test_mat,1)  % iterating over each row in the matrix
    if ~isnan(test_mat(j, 2))
        N = 50; % iterations for animations

        h = 1; gamma = test_mat(j, 1);

        alpha_star = test_mat(j, 2); P_star = test_mat(j, 3);   

        [side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
        P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

        Image = getframe(gcf);
%         folder = 'C:/Users/Work/Desktop/UoA/Year 4/Honours Project/Mathematical Billiards/MATLAB Code/Images/Parallelogram/Animation/Period 2,6,10,14,18,22';
%         baseFileName = sprintf('test%d.jpg', j);
%         fullFileName = fullfile(folder, baseFileName);
%         imwrite(Image.cdata, fullFileName);
    end
end


%% To continue the periodic orbit once it passes through the vertex. We increase the period we want to follow.
% As a period-2 orbit is also a period-6 orbit. We find that at each
% bifurcaton for this periodic orbit, the period increases by 4.

h = 1; gamma = pi/2;    % dimensions of initial parallelogram, we keep h fixed

alpha = pi/2; P = 3.5;  % the angle and position in the square for which we start following from
% (note this normally gives least period-2 orbits in the square)
    
Nmax = 54;

for i=6:4:Nmax
    
    period = i;          % following a period-i orbit
    
    if i == 6   % first iteration
        test_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);
        test_mat = test_mat(1:(end-1),:);   % delete last row
        
        
        %test_mat(:,9) = test_mat(:,4);
        
        % to make all the areas comparable as we have our area func depends on period
        cum_area = test_mat(:,4) - period .* h; % the cumulative area. Area = period*h - cum
        test_mat(:,4) = prod(6:4:Nmax).*h - cum_area .* prod(6:4:Nmax)./period;  
        
        test_mat(:,4) = test_mat(:,4)./(prod(6:4:Nmax).*h)   % normlization, to make numbers smaller
        
        
        
    else % use the previous iterations results to start from
        index = size(test_mat,1);   % last row

        gamma = test_mat(index,1);    % the last gamma that we converged, this is where we start from now
        alpha = test_mat(index,2); P = test_mat(index,3);  % the angle and position in for which we last converged

        temp_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);
        
        %test_mat(:,9) = test_mat(:,4);
        
        % to make all the areas comparable as we have our area func depends on period
        cum_area = temp_mat(:,4) - period .* h; % the cumulative area. Area = period*h - cum
        temp_mat(:,4) = prod(6:4:Nmax).*h - cum_area .* prod(6:4:Nmax)./period;  
        temp_mat(:,4) = temp_mat(:,4)./(prod(6:4:Nmax).*h);   % normlization, to make numbers smaller
        
        test_mat = [test_mat; temp_mat(2:(end-1),:)];  % appending the new results to the matrix

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TEMP PLOTTING

    N = 100; % iterations for animations

    index = size(test_mat,1);
    gamma = test_mat(index, 1);
    alpha_star = test_mat(index, 2);
    P_star = test_mat(index, 3);

    [side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star-0.00, N);
    P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off
    w = waitforbuttonpress;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

%test_mat

%% Saving the matrix

save('Data\parallelogram_period_adding.mat', 'test_mat');

%% For future, we don't need to run the above
test_mat = matfile('Data\parallelogram_period_adding.mat');
test_mat = test_mat.test_mat;


%% Area

plot(test_mat(:,1),test_mat(:,4), 'o', 'Markersize', 3)
xlabel('\gamma', 'FontSize', 22); ylabel('Area', 'FontSize', 22)

set(gca,'XLim',[0 pi/2], 'YLim', [min(test_mat(:,4))-0.02 max(test_mat(:,4))+0.02], 'FontSize', 15)
xticks([0 pi/8 pi/4 3*pi/8 pi/2])
xticklabels({'0', '\pi/8', '\pi/4', '3\pi/8', '\pi/2'})
grid on

% Each sepearte curve is for a different period
% Eg: far most right curve is for period-2, discontinuties are for when we
% hit the vertex

%% Initial Angle

plot(test_mat(:,1),test_mat(:,2), 'o', 'Markersize', 3)
xlabel('\gamma'); ylabel('Initial Angle')

%% Initial Position

plot(test_mat(:,1),test_mat(:,3), 'o', 'Markersize', 4)
xlabel('\gamma'); ylabel('Initial Position')

%% Newton Steps

plot(test_mat(:,1),test_mat(:,5), 'o', 'Markersize', 4)
xlabel('\gamma', 'FontSize', 22); ylabel('Number of Newton Steps', 'FontSize', 22)

set(gca,'XLim',[0 pi/2], 'YLim', [0 max(test_mat(:,5))], 'FontSize', 15)
xticks([0 pi/8 pi/4 3*pi/8 pi/2])
xticklabels({'0', '\pi/8', '\pi/4', '3\pi/8', '\pi/2'})
yticks([0 1 2 3])
yticklabels({'0', '1', '2', '3'})
grid on


%% All the different angles and positions visited 

period = Nmax;

for j = 1:(size(test_mat, 1))    % all entries are good
    gamma = test_mat(j, 1);
    alpha = test_mat(j, 2); P = test_mat(j, 3);
    
    [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
    
    % adding a period column to the matrix test_mat
    next_p = find(abs(p - p(1)) < 1e-6 == 1 & abs(a - a(1)) < 1e-6 == 1, 2); % when do we begin repeating?
    test_mat(j,9) = next_p(2)-next_p(1);
    
    plot3(gamma * ones(size(a)), a, p, 'bo', 'MarkerSize', 1); hold on
    %plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'bo', 'MarkerSize', 1); hold on
end


set(gca,'XLim',[0 pi/2], 'YLim',[0 pi], 'ZLim',[0 8.5], 'FontSize', 15);
xlabel('\gamma', 'FontSize', 22)
ylabel('\alpha_0 - angle', 'FontSize', 18)
zlabel('P_0 - position', 'FontSize', 22)

xticks([0 pi/8 pi/4 3*pi/8 pi/2])
xticklabels({'0', '\pi/8', '\pi/4', '3\pi/8', '\pi/2'})

period = prod(6:4:Nmax);
title(sprintf('Continuation Period-%d', period))

[x y] = meshgrid(0:0.01:pi/2, 0:0.01:pi); 
z0 = ones(size(x));
z1 = 1 + h./sin(x);    
z2 = 2 + h./sin(x);
z3 = 2 + 2*h./sin(x);

surf(x,y,z1,'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on %Plot the surface
surf(x,y,z2,'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on
surf(x,y,z3,'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on
surf(x,y, ones(size(x)),'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on
surf(x,y, zeros(size(x)),'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on

view(0, 0) % az, el

% For dissertation
% x = linspace(0, pi/2, 100);
% plot3(x, zeros(size(x)), zeros(size(x)), 'r-', 'LineWidth', 2); hold on
% plot3(x, zeros(size(x)), ones(size(x)), 'r-', 'LineWidth', 2); hold on
% plot3(x, zeros(size(x)), 1 + h./sin(x), 'r-', 'LineWidth', 2); hold on
% plot3(x, zeros(size(x)), 2 + h./sin(x), 'r-', 'LineWidth', 2); hold on
% plot3(x, zeros(size(x)), 2 + 2*h./sin(x), 'r-', 'LineWidth', 2); hold on




%% Plot of gamma against the period

plot(test_mat(:,1), test_mat(:,9), 'o')
xlabel('\gamma'); ylabel('Period')
xticks([0 pi/4 pi/2]); xticklabels({'0', '\pi/4', '\pi/2'})

set(gca,'XLim',[0 pi/2]);

set(gca,'FontSize', 20)
set(gcf,'color','w'); % background color to white


%% Trying to find empircal formula for the period vs gamma relationship. Note we can use interpolation.

periods = unique(test_mat(:,9));

% When the periodic orbit first and last exists.

last_index = [find(diff(test_mat(:,9))); size(test_mat,1)]; % when do we jump from one period to another?
first_index = [1; last_index(1:(end-1)) + 1];

% Can make a table with this data!
test_mat(first_index,:)

test_mat(last_index,:)

plot(test_mat(first_index,1), test_mat(first_index,9), 'bo'); hold on
plot(test_mat(last_index,1), test_mat(last_index,9), 'ro'); hold on


%%

plot(test_mat(first_index,1), log(test_mat(first_index,9)), 'bo'); % Doesn't seem to exp relationship

%%

plot(test_mat(first_index,1), 1/(test_mat(first_index,9)), 'bo'); % Inverse relationship


%%

plot(log(test_mat(last_index,1)), log(test_mat(last_index,9)), 'bo'); % log log relationship seems to be very good for 
%small gamma


%%

plot(test_mat(first_index,1), log(atan(1./test_mat(first_index,9))), 'bo'); % Inverse relationship


%% TESTING


pp = polyfit(test_mat(last_index,1), log(test_mat(last_index,9)),2);  % Descending powers of polynomials


xq = [0.3:0.01:pi/2];

yq = exp(pp(1) .* xq.^2 + pp(2)* xq + pp(3));

plot(xq, yq, 'go'); hold on
plot(test_mat(last_index,1), test_mat(last_index,9), 'ro')



%%
pp = polyfit(test_mat(last_index,1), test_mat(last_index,9), 5);  % Descending powers of polynomials


xq = [0.3:0.01:0.8]';

yq = xq.^((size(pp, 2)-1):-1:0)   * pp';

plot(xq, yq, 'go'); hold on
plot(test_mat(last_index,1), test_mat(last_index,9), 'ro')


%%

index = 9;


xq = test_mat(last_index(index),1)

disp('TRUE y')
test_mat(last_index(index),9) % TRUE

% APPROXIMATION
disp('APPROX')
xq.^((size(pp, 2)-1):-1:0)   * pp'


%%

A = test_mat(first_index,1).^(0:13);
b = periods;

a = A\b


x = [1];
y = dot(a, x.^(0:13)')

plot(x,y)


%%
xq = [0:0.01:pi/2];

p = polyfit(test_mat(first_index,1), test_mat(first_index,9), 3);

yq = p(1) + p(2) .* xq;


plot(xq, yq, 'o')






%%

%% Investigating the source of the bifurcation

% P vs. g^{6}(alpha, P) plot

h = 1;
P = linspace(0.001, 3.999, 10000);

for i=1:10:size(test_mat,1)
    row_index =  i;

    gamma = test_mat(row_index,1);
    alpha = test_mat(row_index,2);
    period = test_mat(row_index,9);

    for j=1:length(P)
        [side, angle, position] = parallelogram_map(h, gamma, alpha, P(j), period);
        angle_6(j) = angle(period+1); P_6(j) = position(period+1);
    end
    
    plot(P, P_6, 'ro', 'MarkerSize', 2); hold on
    plot(P, P, 'k-'); hold off
    
    xlabel('P', 'FontSize', 22); ylabel('g(\alpha_0, P)', 'FontSize', 22)
    set(gca,'XLim',[0 4], 'YLim',[0 4], 'FontSize', 15);
    title(sprintf('gamma = %0.3f, period = %d', gamma, period))
    Image = getframe(gcf);
end

% Some observations
% Strange droplet behaviour at about gamma = 0.7
% The bifrucations don't occur when we are about to have no more fixed
% points like our period-6 following from earlier. Other points also
% support the periodic orbits, they just don't have maximal area.

%% Trying a surface


h = 1;
row_index =  2000; %i;

gamma = test_mat(row_index,1);
%alpha = test_mat(row_index,2);
period = test_mat(row_index,9);
    
% Our surface

epsilon = 0.001;

[alpha_0 P_0] = meshgrid(0 + epsilon:0.001:pi-epsilon, 0 + epsilon:0.01:1+h/sin(gamma)-epsilon); 
% Generate x = alpha_n and y = P_{n} data
% Not good to include starting angles of exactly 0 and pi as these result
% in discontinuties in our plot.

alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

for row=1:size(alpha_0,1) % over the rows

    for col=1:size(alpha_0,2) % over the columns
        [side_temp alpha_temp pos_temp] = parallelogram_map(h, gamma, alpha_0(row,col), P_0(row,col), period);
        
        ith_iterate = period + 1;

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
surface_no_wall_discont(alpha_0, P_0, P_1, alpha_1);

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

%set(gca,'XLim',[0 4], 'YLim',[0 4], 'FontSize', 15);
title(sprintf('Gamma = %03f, Period = %d', gamma, period))
hold on


% Our cobweb plane, P_n = P_n+period
[x y] = meshgrid(0:0.01:pi, 0:0.1:1+h/sin(gamma)); % Generate x = alpha_n and y = P_{n} data
z = y;     % z = P_{n+period}

s = surf(x,y,z,'FaceAlpha', 0.5, 'FaceColor', 'r'); %Plot the surface
s.EdgeColor = 'none';




%% Trying points, since it's more clear

h = 1;
row_index =  8410; %i;

gamma = test_mat(row_index,1);
%alpha = test_mat(row_index,2);
period = 6; %test_mat(row_index,9);
    
% Our surface

epsilon = 0.001;

[alpha_0 P_0] = meshgrid(0 + epsilon:0.001:pi-epsilon, 0 + epsilon:0.01:1+h/sin(gamma)-epsilon); 
% Generate x = alpha_n and y = P_{n} data
% Not good to include starting angles of exactly 0 and pi as these result
% in discontinuties in our plot.

%alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

counter = 1; % if we find a period-6 orbit

for row=1:size(alpha_0,1) % over the rows

    for col=1:size(alpha_0,2) % over the columns
        [side_temp alpha_temp pos_temp] = parallelogram_map(h, gamma, alpha_0(row,col), P_0(row,col), period);
        
        ith_iterate = period + 1;

        alpha_1(row, col) = alpha_temp(ith_iterate);
        P_1(row, col) = pos_temp(ith_iterate);
        
        % If we found a period-6 point
        if col > 1 & (P_1(row,col-1) < P_0(row,col) < P_1(row,col))
            row
            col
            alpha_0(row,col);
            P_0(row,col);
            break
        end
    end
end





%% Are there any period-6 orbits by using mean value theorem

% If g^(6)(alpha_a, P) - P < 0 and g^(6)(alpha_b, P) - P > 0. Then we
% expect by MVT than for alpha_a < alpha < alpha_b, there is a fixed point.
% This is assuming that g is smooth in this interval.

h = 1; gamma = 0.7;

period = 6;

P_vec = [0 + epsilon:0.01:1+h/sin(gamma)-epsilon];
alpha_vec = [0 + epsilon:0.01:pi-epsilon];

for i=1:length(P_vec)
    P = P_vec(i);
    
    for j=1:length(alpha_vec)
        alpha = alpha_vec(i);
        
        [side_temp alpha_temp pos_temp] = parallelogram_map(h, gamma, alpha_vec(j), P_vec(i), period);
        
        % f^(6)(alpha, P) - alpha
        f(i,j) = alpha_temp(period+1) - alpha;
        
        % g^(6)(alpha, P) - P 
        g(i,j) = pos_temp(period+1) - P;
       
    end
end

%%

index = 189;

x = diff(sign(g(index,:)));

y = find(abs(x) == 2)

% For this value of P
P_vec(index)

sign_change_index = 7;

% between 6th and 7th element, we must cross 0
alpha_vec(y(sign_change_index))
alpha_vec(y(sign_change_index)+1)

%%

P = P_vec(index);   % keep this fixed

alpha_vec2 = linspace(alpha_vec(y(sign_change_index)), alpha_vec(y(sign_change_index)+1), 1000);

for j=1:length(alpha_vec2)
    
    alpha = alpha_vec2(i);
        
    [side_temp alpha_temp pos_temp] = parallelogram_map(h, gamma, alpha_vec2(j), P, period);

    % f^(6)(alpha, P) - alpha
    f_2(j) = alpha_temp(period+1) - alpha;

    % g^(6)(alpha, P) - P 
    g_2(j) = pos_temp(period+1) - P;
    
end

plot(alpha_vec2, f_2, 'o')
figure(2)
plot(alpha_vec2, g_2, 'o')


%%
[a, b] = max(diff(g_2));

alpha_vec2(b+1)

%%
[side, angle, position] = parallelogram_map(h, gamma, alpha_vec2(b+1), P, 9);
P_anim = parallelogram_billiards_draw(h, gamma, side, angle, position);


%%

gamma = 0.7;

[side, angle, position] = parallelogram_map(h, gamma, pi-gamma, 0.5, 4);
P_anim = parallelogram_billiards_draw(h, gamma, side, angle, position);




%% Observing our period-2 to period-6 bifucarion through a surface


h = 1;

for row_index = 1:500:2000  %size(test_mat, 1)
    figure(row_index)
    
    gamma = test_mat(row_index,1);
    period = test_mat(row_index,9);

    % Our surface

    epsilon = 0.001;

    [alpha_0 P_0] = meshgrid(0 + epsilon:0.001:pi-epsilon, 0 + epsilon:0.01:1+h/sin(gamma)-epsilon); 
    % Generate x = alpha_n and y = P_{n} data
    % Not good to include starting angles of exactly 0 and pi as these result
    % in discontinuties in our plot.

    alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

    for row=1:size(alpha_0,1) % over the rows

        for col=1:size(alpha_0,2) % over the columns
            [side_temp alpha_temp pos_temp] = parallelogram_map(h, gamma, alpha_0(row,col), P_0(row,col), period);

            ith_iterate = period + 1;

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
    surface_no_wall_discont(alpha_0, P_0, P_1, alpha_1);

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

    %set(gca,'XLim',[0 4], 'YLim',[0 4], 'FontSize', 15);
    title(sprintf('Gamma = %03f, Period = %d', gamma, period))
    hold on


    % Our cobweb plane, P_n = P_n+period
    [x y] = meshgrid(0:0.01:pi, 0:0.1:1+h/sin(gamma)); % Generate x = alpha_n and y = P_{n} data
    z = y;     % z = P_{n+period}

    s = surf(x,y,z,'FaceAlpha', 0.5, 'FaceColor', 'r'); hold off; %Plot the surface
    s.EdgeColor = 'none';
    Image = getframe(gcf);
end



%% Period-2 to period-6 up to... period-54 saved data.

test_mat = matfile('Data\parallelogram_period_adding.mat');
test_mat = test_mat.test_mat;

%% Period subtraction. Period-54 to period-50, period ..., running the parameter gamma backwards

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make sure that the 'gamma_perturb' argument in
% parallelogram_newton_solver.m is negative in order to reverse gamma
% direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note that we collapse and all the points in the trajectory collapse into
% a vertex! This can be thought of as a period-6 orbit I guess.

h = 1; gamma = test_mat(end,1);    % dimensions of initial parallelogram, we keep h fixed

alpha = test_mat(end,2); P = test_mat(end,3);  % the angle and position in the square for which we start following from
    
Nmax = 54;

for i=Nmax:-4:10
    
    period = i;          % following a period-i orbit
    
    if i == Nmax   % first iteration
        test_mat2 = parallelogram_newton_solver(h, gamma, alpha, P, period);
        test_mat2 = test_mat2(1:(end-1),:);   % delete last row
        
        % to make all the areas comparable as we have our area func depends on period
        cum_area = test_mat2(:,4) - period .* h; % the cumulative area. Area = period*h - cum
        test_mat2(:,4) = prod(6:4:Nmax).*h - cum_area .* prod(6:4:Nmax)./period;  
        
        test_mat2(:,4) = test_mat2(:,4)./(prod(6:4:Nmax).*h)   % normlization, to make numbers smaller
        
        test_mat2(:,9) = period*ones(size(test_mat2,1),1); % the period
        
    else % use the previous iterations results to start from
        index = size(test_mat2,1);   % last row

        gamma = test_mat2(index,1);    % the last gamma that we converged, this is where we start from now
        alpha = test_mat2(index,2); P = test_mat2(index,3);  % the angle and position in for which we last converged

        temp_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);
        
        %test_mat(:,9) = test_mat(:,4);
        
        % to make all the areas comparable as we have our area func depends on period
        cum_area = temp_mat(:,4) - period .* h; % the cumulative area. Area = period*h - cum
        temp_mat(:,4) = prod(6:4:Nmax).*h - cum_area .* prod(6:4:Nmax)./period;  
        temp_mat(:,4) = temp_mat(:,4)./(prod(6:4:Nmax).*h);   % normlization, to make numbers smaller
        
        temp_mat(:,9) = period*ones(size(temp_mat,1),1); % the period
        
        test_mat2 = [test_mat2; temp_mat(2:(end-1),:)];  % appending the new results to the matrix
    
    temp_mat
    period
    max(temp_mat(1:(end-1),5))
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TEMP PLOTTING

    N = 100; % iterations for animations

    index = size(test_mat2,1);
    gamma = test_mat2(index, 1);
    alpha_star = test_mat2(index, 2);
    P_star = test_mat2(index, 3);

    [side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star-0.00, N);
    P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off
    w = waitforbuttonpress;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

%%

save('Data\parallelogram_period_subtracting.mat', 'test_mat2');

%% For future, we don't need to run the above
test_mat2 = matfile('Data\parallelogram_period_subtracting.mat');
test_mat2 = test_mat2.test_mat2;

%% Animating above

for j=1:20:size(test_mat2,1)  % iterating over each row in the matrix
    if ~isnan(test_mat2(j, 2))
        N = 50; % iterations for animations

        h = 1; gamma = test_mat2(j, 1);

        alpha_star = test_mat2(j, 2); P_star = test_mat2(j, 3);   

        [side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
        P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

        Image = getframe(gcf);
%         folder = 'C:/Users/Work/Desktop/UoA/Year 4/Honours Project/Mathematical Billiards/MATLAB Code/Images/Parallelogram/Animation/Period 54,50,46,42,46';
%         baseFileName = sprintf('test%d.jpg', j);
%         fullFileName = fullfile(folder, baseFileName);
%         imwrite(Image.cdata, fullFileName);
    end
end

%%
period = Nmax;

for j = 1:(size(test_mat2, 1))    % all entries are good
    gamma = test_mat2(j, 1);
    alpha = test_mat2(j, 2); P = test_mat2(j, 3);
    
    [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
    
    plot3(gamma * ones(size(a)), a, p, 'bo', 'MarkerSize', 1); hold on
    %plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'bo', 'MarkerSize', 1); hold on
end


set(gca,'XLim',[0 pi/2], 'YLim',[0 pi], 'ZLim',[0 8.5], 'FontSize', 15);
xlabel('\gamma', 'FontSize', 22)
ylabel('\alpha_0 - angle', 'FontSize', 18)
zlabel('P_0 - position', 'FontSize', 22)

xticks([0 pi/8 pi/4 3*pi/8 pi/2])
xticklabels({'0', '\pi/8', '\pi/4', '3\pi/8', '\pi/2'})

period = prod(10:4:Nmax);
title(sprintf('Continuation Period-%d', period))

[x y] = meshgrid(0:0.01:pi/2, 0:0.01:pi); 
z0 = ones(size(x));
z1 = 1 + h./sin(x);    
z2 = 2 + h./sin(x);
z3 = 2 + 2*h./sin(x);

surf(x,y,z1,'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on %Plot the surface
surf(x,y,z2,'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on
surf(x,y,z3,'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on
surf(x,y, ones(size(x)),'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on
surf(x,y, zeros(size(x)),'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on

view(0, 0) % az, el

% For dissertation
% x = linspace(0, pi/2, 100);
% plot3(x, zeros(size(x)), zeros(size(x)), 'r-', 'LineWidth', 2); hold on
% plot3(x, zeros(size(x)), ones(size(x)), 'r-', 'LineWidth', 2); hold on
% plot3(x, zeros(size(x)), 1 + h./sin(x), 'r-', 'LineWidth', 2); hold on
% plot3(x, zeros(size(x)), 2 + h./sin(x), 'r-', 'LineWidth', 2); hold on
% plot3(x, zeros(size(x)), 2 + 2*h./sin(x), 'r-', 'LineWidth', 2); hold on





%% Comparing the forward and backwards for co-existence of periodic orbits

test_mat = matfile('Data\parallelogram_period_adding.mat');
test_mat = test_mat.test_mat;

test_mat2 = matfile('Data\parallelogram_period_subtracting.mat');
test_mat2 = test_mat2.test_mat2;

Nmax = 54;
h = 1;

period = Nmax;

for j = 1:10:(size(test_mat, 1))    % all entries are good
    gamma = test_mat(j, 1);
    alpha = test_mat(j, 2); P = test_mat(j, 3);
    period = test_mat(j, 9);
    
    [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
    
    %plot3(gamma * ones(size(a)), a, p, 'bo', 'MarkerSize', 1); hold on
    %plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'bo', 'MarkerSize', 1); hold on
    
    plot(gamma * ones(size(a)), p, 'bo', 'MarkerSize', 1); hold on % Just the gamma vs. position
end

waitforbuttonpress
for j = 1:50:(size(test_mat2, 1))    % all entries are good
    gamma = test_mat2(j, 1);
    alpha = test_mat2(j, 2); P = test_mat2(j, 3);
    period = test_mat2(j, 9);
    
    [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
    
    %plot3(gamma * ones(size(a)), a, p, 'bo', 'MarkerSize', 1); hold on
    %plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'bo', 'MarkerSize', 1); hold on
    
    plot(gamma * ones(size(p)), p, 'r*', 'MarkerSize', 1); hold on % Just the gamma vs. position
end
waitforbuttonpress
set(gca,'XLim',[0 pi/2], 'YLim',[0 8.5], 'FontSize', 15);
xlabel('\gamma', 'FontSize', 22)
ylabel('P_0 - position', 'FontSize', 22)

xticks([0 pi/8 pi/4 3*pi/8 pi/2])
xticklabels({'0', '\pi/8', '\pi/4', '3\pi/8', '\pi/2'})

period = prod(10:4:Nmax);
title(sprintf('Continuation Period-%d', period))

% For dissertation
x = linspace(0, pi/2, 100);
plot(x, zeros(size(x)), 'k-', 'LineWidth', 2); hold on
plot(x, ones(size(x)), 'k-', 'LineWidth', 2); hold on
plot(x, 1 + h./sin(x), 'k-', 'LineWidth', 2); hold on
plot(x, 2 + h./sin(x), 'k-', 'LineWidth', 2); hold on
plot(x, 2 + 2*h./sin(x), 'k-', 'LineWidth', 2); hold on







%% Are there any least period-4 orbits in the parallelogram?

% I don't think so based on the angle not returning as the same. Also since
% our continuation just fails straight away.

% Note that we have the least period-2 intersection here

h = 1;

gamma = 0.7; % random angle
period = 4;
    
% Our surface

epsilon = 0.001;

[alpha_0 P_0] = meshgrid(0 + epsilon:0.001:pi-epsilon, 0 + epsilon:0.01:1+h/sin(gamma)-epsilon); 
% Generate x = alpha_n and y = P_{n} data
% Not good to include starting angles of exactly 0 and pi as these result
% in discontinuties in our plot.

alpha_1 = zeros(size(P_0));   P_1 = zeros(size(P_0));

for row=1:size(alpha_0,1) % over the rows

    for col=1:size(alpha_0,2) % over the columns
        [side_temp alpha_temp pos_temp] = parallelogram_map(h, gamma, alpha_0(row,col), P_0(row,col), period);
        
        ith_iterate = period + 1;

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

%set(gca,'XLim',[0 4], 'YLim',[0 4], 'FontSize', 15);
title(sprintf('Gamma = %03f, Period = %d', gamma, period))
hold on


% Our cobweb plane, P_n = P_n+period
[x y] = meshgrid(0:0.01:pi, 0:0.1:1+h/sin(gamma)); % Generate x = alpha_n and y = P_{n} data
z = y;     % z = P_{n+period}

s = surf(x,y,z,'FaceAlpha', 0.5, 'FaceColor', 'r') %Plot the surface
s.EdgeColor = 'none';





%% NEW JACOBIAN TESTING. (solving the alpha-P equation (singular))
alpha_star = 1.324; P_star = 0.123; period = 4;

h = 1.236; gamma = pi/2-0.1;

parallelogram_jacobian(h, gamma, alpha_star, P_star, period)

% ONLY ONE EIGENVALUE 1 now (compared to square and rectangle)

%%
clear all; clc;

h = 1; gamma = pi/2-0.1;


init_alpha = pi/3; init_P = 0.1;
init_guess = [init_alpha; init_P]; % init_angle and init_pos

period = 4; % what period are we seaching for?

alpha = init_alpha; P = init_P;
guess=  init_guess;

for j=1:7
    disp(j)
    % F^{N}(\alpha, \P) term
   [F_alpha, F_P] = parallelogram_map(h, gamma, alpha, P, period);
    F_n = [F_alpha(period+1); F_P(period+1)];
    
    parallelogram_jacobian(h, gamma, alpha, P, period)- eye(2)
    
   guess = guess - inv(parallelogram_jacobian(h, gamma, alpha, P, period) - eye(2)) * (F_n - guess);
   
   alpha = mod(guess(1), pi); P = mod(guess(2), 1);
   
end