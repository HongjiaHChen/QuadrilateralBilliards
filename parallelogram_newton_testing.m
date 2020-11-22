%% All the different initial positions and angles which give a period-6 orbit in the square and maximise the area
% result in the same continuation!

h = 1;
gamma = pi/2;
period = 6;

init_alpha_1 = atan(2);

bigmat_0p25 = parallelogram_newton_solver(h, gamma, init_alpha_1, 0.25, period);
bigmat_0p75 = parallelogram_newton_solver(h, gamma, init_alpha_1, 0.75, period);
bigmat_2p25 = parallelogram_newton_solver(h, gamma, init_alpha_1, 2.25, period);
bigmat_2p75 = parallelogram_newton_solver(h, gamma, init_alpha_1, 2.75, period);

bigmat_1p25 = parallelogram_newton_solver(h, gamma, init_alpha_1, 1.25, period); % a different P.O in the square
bigmat_1p75 = parallelogram_newton_solver(h, gamma, init_alpha_1, 1.75, period); % a different P.O in the square

init_alpha_2 = atan(1/2);

bigmat_1p5 = parallelogram_newton_solver(h, gamma, init_alpha_2, 1.5, period);
bigmat_3p5 = parallelogram_newton_solver(h, gamma, init_alpha_2, 3.5, period);


%%

for j = 1:(size(bigmat_0p25, 1)-1)    % last entry is the one we terminate at, so is scuffed
    gamma = bigmat_0p25(j, 1);
    alpha = bigmat_0p25(j, 2); P = bigmat_0p25(j, 3);
    
    [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
    
    plot3(gamma * ones(size(a)), a, p, 'bo', 'MarkerSize', 1); hold on
    %plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'bo', 'MarkerSize', 1); hold on
end
view(0,0)
w = waitforbuttonpress; 

for j = 1:(size(bigmat_1p25, 1)-1)    % last entry is the one we terminate at, so is scuffed
    gamma = bigmat_1p25(j, 1);
    alpha = bigmat_1p25(j, 2); P = bigmat_1p25(j, 3);
    
    [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
    
    plot3(gamma * ones(size(a)), a, p, 'bo', 'MarkerSize', 1); hold on
    %plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'bo', 'MarkerSize', 1); hold on
end
view(0,0)
w = waitforbuttonpress; 

for j = 1:(size(bigmat_0p75, 1)-1)    % last entry is the one we terminate at, so is scuffed
    gamma = bigmat_0p75(j, 1);
    alpha = bigmat_0p75(j, 2); P = bigmat_0p75(j, 3);
    
    [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
    
    plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'go', 'MarkerSize', 2); hold on
end
view(0,0)
w = waitforbuttonpress; 

for j = 1:(size(bigmat_2p25, 1)-1)    % last entry is the one we terminate at, so is scuffed
    gamma = bigmat_2p25(j, 1);
    alpha = bigmat_2p25(j, 2); P = bigmat_2p25(j, 3);
    
    [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
    
    plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'c*', 'MarkerSize', 1); hold on
end
view(0,0)
w = waitforbuttonpress; 

for j = 1:(size(bigmat_2p75, 1)-1)    % last entry is the one we terminate at, so is scuffed
    gamma = bigmat_2p75(j, 1);
    alpha = bigmat_2p75(j, 2); P = bigmat_2p75(j, 3);
    
    [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
    
    plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'm+', 'MarkerSize', 1); hold on
end

for j = 1:(size(bigmat_3p5, 1)-1)    % last entry is the one we terminate at, so is scuffed
    gamma = bigmat_3p5(j, 1);
    alpha = bigmat_3p5(j, 2); P = bigmat_3p5(j, 3);
    
    [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
    
    plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'ro', 'MarkerSize', 1); hold on
end



for j = 1:(size(bigmat_1p25, 1)-1)    % last entry is the one we terminate at, so is scuffed
    gamma = bigmat_1p25(j, 1);
    alpha = bigmat_1p25(j, 2); P = bigmat_1p25(j, 3);
    
    [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
    
    plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'ko', 'MarkerSize', 1); hold on
end
% 
% 
% for j = 1:(size(bigmat_1p75, 1)-1)    % last entry is the one we terminate at, so is scuffed
%     gamma = bigmat_1p75(j, 1);
%     alpha = bigmat_1p75(j, 2); P = bigmat_1p75(j, 3);
%     
%     [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
%     
%     plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'yo', 'MarkerSize', 1); hold on
% end

set(gca,'XLim',[1.35 pi/2], 'YLim',[0 pi], 'ZLim',[0 4]);
xlabel('Gamma of parallelogram')
ylabel('\alpha_0 - angle')
zlabel('P_0 - position')
title(sprintf('Continuation Period-%d', period))

view(0,0)

%% Clearly it appears that the positions visited are exactly the same and we terminate at the same point.
% In fact, P_0 = 0.25 and P_0 = 2.75 in the square are the exact same
% trajectory
% whereas  P_0 = 0.75 and P_0 = 2.25 are the exact same trajectory. It is
% in fact the same as above except in the opposite direction.
% The angle differs by pi-angle for these two groups

% view(90,0) % to see the above

% What do the orbits in the billiard look like? % They look the same!
N = 10; % iterations for animations


index = 1800 % what row of the matrix do we want?

% For P_0 = 0.25,
gamma = bigmat_0p25(index, 1);
alpha_star = bigmat_0p25(index, 2);
P_star = bigmat_0p25(index, 3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

w = waitforbuttonpress;  

% For P_0 = 2.75,
gamma = bigmat_2p75(index, 1);
alpha_star = bigmat_2p75(index, 2);
P_star = bigmat_2p75(index, 3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

w = waitforbuttonpress;  

% For P_0 = 0.75,
gamma = bigmat_0p75(index, 1);
alpha_star = bigmat_0p75(index, 2);
P_star = bigmat_0p75(index, 3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

w = waitforbuttonpress;  





%% Is there a saddle node bifurcation? Why does the orbit dissapear?

% Forward direction
h = 1;
gamma = pi/2;
period = 6;

init_alpha_1 = atan(2);

bigmat_0p75 = parallelogram_newton_solver(h, gamma, init_alpha_1, 0.75, period);

%% Backwards now
index = size(bigmat_0p75,1)-1; % last good row

% For P_0 = 0.75,
gamma = bigmat_0p75(index, 1);
alpha_star = bigmat_0p75(index, 2);
P_star = bigmat_0p75(index, 3)+0.00; % bump this up a bit

bigmat_0p75_back = parallelogram_newton_solver(h, gamma, alpha_star, P_star, period);

%% Plot of the positions and angles vs gamma

for j = 1:(size(bigmat_0p75_back, 1)-1)    % last entry is the one we terminate at, so is scuffed
    gamma = bigmat_0p75_back(j, 1);
    alpha = bigmat_0p75_back(j, 2); P = bigmat_0p75_back(j, 3);
    
    [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
    
    plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'go', 'MarkerSize', 2); hold on
end

xlabel('\gamma'); ylabel('\alpha_0 - angle');
zlabel('P_0 - position')



%% What values of P support a period-6 orbit?
h = 1;
period = 6;

row_index =  1000;
bigmat_0p25(row_index,:)

% Fix gamma
gamma = bigmat_0p25(row_index,1);
% Fix alpha_0
alpha = bigmat_0p25(row_index,2);

%P = linspace(0.001, 0.999, 10000);
P = linspace(0.001, 3.999, 10000);

for j=1:length(P)
    
    [side, angle, position] = parallelogram_map(h, gamma, alpha, P(j), period);
    angle_6(j) = angle(period+1); P_6(j) = position(period+1);
    
end

%%

plot(P, P_6, 'o', 'MarkerSize', 2); hold on
plot(P, P, 'r-');
xlabel('P', 'FontSize', 22); ylabel('g(\alpha_0, P)', 'FontSize', 22)
set(gca,'XLim',[0 1], 'YLim',[0 1], 'FontSize', 15);


%% Looping over the different gammas

for i=1:100:size(bigmat_0p25,1)
    row_index =  i;

    gamma = bigmat_0p25(row_index,1);
    alpha = bigmat_0p25(row_index,2);

    for j=1:length(P)
        [side, angle, position] = parallelogram_map(h, gamma, alpha, P(j), period);
        angle_6(j) = angle(period+1); P_6(j) = position(period+1);
    end
    
    plot(P, P_6, 'ro', 'MarkerSize', 2); hold on
    plot(P, P, 'k-'); hold off
    
    xlabel('P', 'FontSize', 22); ylabel('g(\alpha_0, P)', 'FontSize', 22)
    set(gca,'XLim',[0 4], 'YLim',[0 4], 'FontSize', 15);
    title(sprintf('gamma = %0.3f', gamma))
    Image = getframe(gcf);
end

% Looks like there will be another Period-6 point soon, although it will
% simply pass through


%% How does this look like for P_0=1.25

for i=1:10:size(bigmat_1p25,1)
    row_index =  i;

    gamma = bigmat_1p25(row_index,1);
    alpha = bigmat_1p25(row_index,2);

    for j=1:length(P)
        [side, angle, position] = parallelogram_map(h, gamma, alpha, P(j), period);
        angle_6(j) = angle(period+1); P_6(j) = position(period+1);
    end
    
    plot(P, P_6, 'ro', 'MarkerSize', 2); hold on
    plot(P, P, 'k-'); hold off
    
    xlabel('P', 'FontSize', 22); ylabel('g(\alpha_0, P)', 'FontSize', 22)
    set(gca,'XLim',[0 4], 'YLim',[0 4], 'FontSize', 15);
    title(sprintf('gamma = %0.3f', gamma))
    Image = getframe(gcf);
end

% This one appears different to the other one, we have the red map passing
% through the y=x as well as an interesting droplet formation




%% A surface

h = 1;
row_index =  2112; %i;

gamma = bigmat_1p25(row_index,1);
%alpha = test_mat(row_index,2);
period = 6;
    
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

        ith_iterate = 4;   %depends on the N we have chosen , x_1, x_2, x_3,...      
        % change this for any value from 2 to N+1, 
        % for a period  of ith_iterate - 1
        
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

%set(gca,'XLim',[0 pi], 'YLim',[0 1+h/sin(gamma)], 'ZLim',[0 1+h/sin(gamma)], 'FontSize', 15);
set(gca,'XLim',[1.2 1.4], 'YLim',[0 1+h/sin(gamma)], 'ZLim',[0 1+h/sin(gamma)], 'FontSize', 15);
title(sprintf('Gamma = %03f, Period = %d', gamma, period))
hold on


% Our cobweb plane, P_n = P_n+period
[x y] = meshgrid(0:0.01:pi, 0:0.1:1+h/sin(gamma)); % Generate x = alpha_n and y = P_{n} data
z = y;     % z = P_{n+period}

s = surf(x,y,z,'FaceAlpha', 0.5, 'FaceColor', 'r'); %Plot the surface
s.EdgeColor = 'none';



%%

gamma = test_mat2(end-200,1); %test_mat(end,1);  %test_mat2(end-50,1);
alpha = test_mat2(end-200,2);%test_mat(end,2); %test_mat2(end-50,2);
P = test_mat2(end-200,3);    %test_mat(end,3);%test_mat2(end-50,3);

[side, angle, position] = parallelogram_map(h, gamma, alpha, P, 100);
P_anim = parallelogram_billiards_draw(h, gamma, side, angle, position);












%% TESTING DELETE











%% TRYING TO FOLLOW SOME MORE

j = (size(bigmat_0p25, 1)-1);
gamma = bigmat_0p25(j, 1);
alpha = bigmat_0p25(j, 2); P = bigmat_0p25(j, 3);
[s, a, p] = parallelogram_map(h, gamma, alpha, P, period);


index = 5;

test_mat = parallelogram_newton_solver(h, gamma, a(index)+0.01, p(index)+0.01, period);

test_mat




%% FOLLOWING REPEATED Period-6
row = 850;
bigmat(row,:)

%% Testing, feel free to change this
N = 8;

h = 1; gamma = bigmat(row, 1) ; 
alpha = bigmat(row, 2); P = bigmat(row, 3);  

test_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);

test_mat



%%

index = 200; % size(test_mat,1)-10;
gamma = test_mat(index, 1);
alpha_star = test_mat(index, 2);
P_star = test_mat(index, 3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);















%% Following randoms

period = 96;
h = 1; gamma = pi/2; %bigmat(row, 1) ; 
alpha = pi/2; P = 3.5;  

test_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);

test_mat

%%
N = 100; % iterations for animations

index = size(test_mat,1)-1;
gamma = test_mat(index, 1);
alpha_star = test_mat(index, 2);
P_star = test_mat(index, 3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star-0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);

%%
index = size(test_mat,1)-1

alpha_star = test_mat(index, 2); P_star = test_mat(index, 3)-0.4;

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);




%%
























%% TESTING OTHER PERIOD-2

clear all; clc
tic
period = 6*1; % what period are we seaching for?

h = 1; % keep this fixed for now

% We know for square that the angle and position that is a period-6 orbit
% and that maximises the area is angle=atan(2), pos = 0.25. This is our
% initial guess.

gamma = pi/2;             % for the square
alpha = pi/2; P = 2.5;   % Change the P value

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







%% ANOTHER PERIOD-6

clear all; clc
period = 2*12; % what period are we seaching for?

h = 1; % keep this fixed for now

gamma = pi/2;             % for the square
alpha = pi/2; P = 1.5;

test_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);
test_mat = test_mat(1:(end-1),:)



%%

index = size(test_mat,1);   % last row

gamma = test_mat(index,1);    % the last gamma that we converged, this is where we start from now
alpha = test_mat(index,2); P = test_mat(index,3);  % the angle and position in for which we last converged

period = 18;

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

row = size(test_mat,1) - 10;
h = 1; gamma = test_mat(row, 1) ;

alpha_star = test_mat(row, 2); P_star = test_mat(row, 3);   
[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star+0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);



%% Animating the continuation of orbits

for j=600:10:size(test_mat,1)  % iterating over each row in the matrix
    if ~isnan(test_mat(j, 2))
        N = 50; % iterations for animations

        h = 1; gamma = test_mat(j, 1);

        alpha_star = test_mat(j, 2); P_star = test_mat(j, 3);   

        [side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
        
        % adding a period column to the matrix test_mat
        next_p = find(abs(position - position(1)) < 1e-6 == 1 & abs(alpha - alpha(1)) < 1e-6 == 1, 2); % when do we begin repeating?
        test_mat(j,9) = next_p(2)-next_p(1);
    
        P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

        Image = getframe(gcf);
%         folder = 'C:/Users/Work/Desktop/UoA/Year 4/Honours Project/Mathematical Billiards/MATLAB Code/Images/Parallelogram/Animation/Period 6';
%         baseFileName = sprintf('test%d.jpg', j);
%         fullFileName = fullfile(folder, baseFileName);
%         imwrite(Image.cdata, fullFileName);
    end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Following larger periodic orbits

%% Period 8. alpha_0 = atan(3), P_0 = 0.5
h = 1;
gamma = pi/2;
period = 8;

init_alpha_1 = atan(3); init_P_1 = 0.5;

p8_alpha_1 = parallelogram_newton_solver(h, gamma, init_alpha_1, init_P_1, period);
p8_alpha_1(end,:)

%% Period 8. alpha_0 = atan(1/3), P_0 = 0.5
h = 1;
gamma = pi/2;
period = 8;

init_alpha_2 = atan(1/3); init_P_2 = 0.5;

p8_alpha_2 = parallelogram_newton_solver(h, gamma, init_alpha_2, init_P_2, period);
p8_alpha_2(end,:)

%% Period 10. alpha_0 = atan(4), P_0 = 1/8
h = 1;
gamma = pi/2;
period = 10;

init_alpha_1 = atan(4); init_P_1 = 1/8;

p8_alpha_1 = parallelogram_newton_solver(h, gamma, init_alpha_1, init_P_1, period);
p8_alpha_1(end,:)

%% Period 10. alpha_0 = atan(1/4), P_0 = 1/2
h = 1;
gamma = pi/2;
period = 10;

init_alpha_1 = atan(1/4); init_P_1 = 1/2;

p8_alpha_1 = parallelogram_newton_solver(h, gamma, init_alpha_1, init_P_1, period);
p8_alpha_1(end,:)

%% Period 10. alpha_0 = atan(3/2), P_0 = 1/2
h = 1;
gamma = pi/2;
period = 10;

init_alpha_1 = atan(3/2); init_P_1 = 1/2;

p8_alpha_1 = parallelogram_newton_solver(h, gamma, init_alpha_1, init_P_1, period);
p8_alpha_1(end,:)

%% Period 10. alpha_0 = atan(2/3), P_0 = 1/4
h = 1;
gamma = pi/2;
period = 10;

init_alpha_1 = atan(2/3); init_P_1 = 1/4;

p8_alpha_1 = parallelogram_newton_solver(h, gamma, init_alpha_1, init_P_1, period);
p8_alpha_1(end,:)


