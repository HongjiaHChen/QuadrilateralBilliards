% Dissertation images for square to parallelogram continuation



%% Pure period-6 continution than results in termination after vertex intersection
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


%% Animation

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

%% Choosing four specific frames for dissertation

% We want gamma = pi/2, gamma = ..., gamma = ..., gamma = 1.3651
% The inbetween gammas, we pick to be 1/3 of the way between 1.3651 and
% pi/2 and 2/3s of the way.
N = 10;

%% For gamma  = pi/2,

gamma = test_mat(1, 1);
alpha_star = test_mat(1, 2); P_star = test_mat(1,3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

text(-0.3, 1, '(a)', 'Color','k', 'FontSize', 22)

%% For gamma = ..., the 6858th row in the matrix. Found by computing size(test_mat, 1)/3

row_index = round(size(test_mat, 1)/3, 0)

gamma = test_mat(row_index, 1)
alpha_star = test_mat(row_index, 2); P_star = test_mat(row_index,3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

text(-0.3, 1, '(c)', 'Color','k', 'FontSize', 22)

%% For gamma = ..., the 13716th row in the matrix. Found by computing size(test_mat, 1)/3

row_index = round(2*size(test_mat, 1)/3, 0)

gamma = test_mat(row_index, 1)
alpha_star = test_mat(row_index, 2); P_star = test_mat(row_index,3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

text(-0.3, 1, '(e)', 'Color','k', 'FontSize', 22)


%% For gamma = 1.3651

gamma = test_mat(end, 1)
alpha_star = test_mat(end, 2); P_star = test_mat(end,3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

text(-0.3, 1, '(g)', 'Color','k', 'FontSize', 22)

%% Which positions support these period-6 points?

bigmat_0p25 = test_mat;

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

%% Looping over the different gammas

bigmat_0p25 = test_mat;

for i=1:100:size(bigmat_0p25,1)-1
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


%% The apparent intersection for 1<P<2 is at: P = 1.15, when alpha = 1.3129
% The P matches the same after 6 iterations, however the angle is
% different.

gamma = test_mat(end, 1);
alpha_star = test_mat(end, 2); P_star = 1.15;

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, 7);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off



%% For disseration, we only need the four specific gammas

%% For gamma = pi/2,

row_index = 1;

gamma = bigmat_0p25(row_index,1);
alpha = bigmat_0p25(row_index,2);

for j=1:length(P)
    [side, angle, position] = parallelogram_map(h, gamma, alpha, P(j), period);
    angle_6(j) = angle(period+1); P_6(j) = position(period+1);
end

set(gca,'fontsize',140, 'FontWeight', 'Bold')

plot(P, P, 'b-', 'LineWidth', 4); pbaspect([1 1 1]); hold on
plot(P, P_6, 'ro', 'MarkerSize', 4)

xlabel('$P_0$','Interpreter','latex','String','{\boldmath$P_0$}'); 
ylabel('$g^{(6)}(\alpha_0, P_0)$', 'Interpreter','latex','String', '{\boldmath$g^{(6)}(\alpha_0, P_0)$}')
xticks([0 1 1+1/sin(gamma)]); xticklabels({'A', 'B', 'C'})
yticks([0 1 1+1/sin(gamma)]); yticklabels({'A', 'B', 'C'})

set(gca,'XLim',[0 1+1/sin(gamma)], 'YLim',[0 1+1/sin(gamma)], 'FontSize', 15);
%title(sprintf('gamma = %0.3f', gamma))
set(gcf,'color','w');

text(-0.3, 1+1/sin(gamma), '(b)', 'Color','k', 'FontSize', 22)


%% For gamma = ..., 1/3 of the way

row_index = round(size(test_mat, 1)/3,0);

gamma = bigmat_0p25(row_index,1);
alpha = bigmat_0p25(row_index,2);

for j=1:length(P)
    [side, angle, position] = parallelogram_map(h, gamma, alpha, P(j), period);
    angle_6(j) = angle(period+1); P_6(j) = position(period+1);
end

set(gca,'fontsize',140, 'FontWeight', 'Bold')

plot(P, P, 'b-', 'LineWidth', 4); pbaspect([1 1 1]); hold on
plot(P, P_6, 'ro', 'MarkerSize', 4)

xlabel('$P_0$','Interpreter','latex','String','{\boldmath$P_0$}', 'Fontsize', 40); 
ylabel('$g^{(6)}(\alpha_0, P_0)$', 'Interpreter','latex','String', '{\boldmath$g^{(6)}(\alpha_0, P_0)$}')
xticks([0 1 1+1/sin(gamma)]); xticklabels({'A', 'B', 'C'})
yticks([0 1 1+1/sin(gamma)]); yticklabels({'A', 'B', 'C'})

set(gca,'XLim',[0 1+1/sin(gamma)], 'YLim',[0 1+1/sin(gamma)], 'FontSize', 15);
%title(sprintf('gamma = %0.3f', gamma))
set(gcf,'color','w');

text(-0.3, 1+1/sin(gamma), '(d)', 'Color','k', 'FontSize', 22)


%% For gamma = ..., 2/3 of the way

row_index = round(2*size(test_mat, 1)/3,0);

gamma = bigmat_0p25(row_index,1);
alpha = bigmat_0p25(row_index,2);

for j=1:length(P)
    [side, angle, position] = parallelogram_map(h, gamma, alpha, P(j), period);
    angle_6(j) = angle(period+1); P_6(j) = position(period+1);
end

set(gca,'fontsize',140, 'FontWeight', 'Bold')

plot(P, P, 'b-', 'LineWidth', 4); hold on 
plot(P, P_6, 'ro', 'MarkerSize', 4); pbaspect([1 1 1])

set(gca,'fontsize',140, 'FontWeight', 'Bold')

xlabel('$P_0$','Interpreter','latex','String','{\boldmath$P_0$}', 'Fontsize', 40); 
ylabel('$g^{(6)}(\alpha_0, P_0)$', 'Interpreter','latex','String', '{\boldmath$g^{(6)}(\alpha_0, P_0)$}')
xticks([0 1 1+1/sin(gamma)]); xticklabels({'A', 'B', 'C'})
yticks([0 1 1+1/sin(gamma)]); yticklabels({'A', 'B', 'C'})
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a, 'fontsize',100)

set(gca,'XLim',[0 1+1/sin(gamma)], 'YLim',[0 1+1/sin(gamma)], 'FontSize', 15);
%title(sprintf('gamma = %0.3f', gamma))
set(gcf,'color','w');

text(-0.3, 1+1/sin(gamma), '(f)', 'Color','k', 'FontSize', 22)

%% For gamma = 1.3651

row_index = size(test_mat,1)-1;

gamma = bigmat_0p25(row_index,1);
alpha = bigmat_0p25(row_index,2);

for j=1:length(P)
    [side, angle, position] = parallelogram_map(h, gamma, alpha, P(j), period);
    angle_6(j) = angle(period+1); P_6(j) = position(period+1);
end

set(gca,'fontsize',140, 'FontWeight', 'Bold')

plot(P, P, 'b-', 'LineWidth', 4); pbaspect([1 1 1]); hold on
plot(P, P_6, 'ro', 'MarkerSize', 4)

xlabel('$P_0$','Interpreter','latex','String','{\boldmath$P_0$}'); 
ylabel('$g^{(6)}(\alpha_0, P_0)$', 'Interpreter','latex','String', '{\boldmath$g^{(6)}(\alpha_0, P_0)$}')
xticks([0 1 1+1/sin(gamma)]); xticklabels({'A', 'B', 'C'})
yticks([0 1 1+1/sin(gamma)]); yticklabels({'A', 'B', 'C'})

set(gca,'XLim',[0 1+1/sin(gamma)], 'YLim',[0 1+1/sin(gamma)], 'FontSize', 15);
%title(sprintf('gamma = %0.3f', gamma))
set(gcf,'color','w');

text(-0.3, 1+1/sin(gamma), '(h)', 'Color','k', 'FontSize', 22)



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Period adding bifucation figures

h = 1;

% Loading the saved data
test_mat = matfile('Data\parallelogram_period_adding.mat');
test_mat = test_mat.test_mat;


%% Static image of period adding bifurcations

N = 100; % iterations for animations

%% Initially

index = 1;
gamma = test_mat(index, 1);
alpha_star = test_mat(index, 2);
P_star = test_mat(index, 3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star-0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off
    
text(-0.25, 1, '(a)', 'FontSize', 22)

%% Period-2 about to hit vertex

index = find(test_mat(:,9) == 2, 1, 'last')-1000;
gamma = test_mat(index, 1);
alpha_star = test_mat(index, 2);
P_star = test_mat(index, 3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star-0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off
    
text(0, 1, '(b)', 'FontSize', 22)


%% Formation of Period-6 orbit

index = find(test_mat(:,9) == 6, 1, 'first');
gamma = test_mat(index, 1);
alpha_star = test_mat(index, 2);
P_star = test_mat(index, 3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star-0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off
    
text(0, 1, '(c)', 'FontSize', 22)

%% Period-6 orbit about to hit vertex

index = find(test_mat(:,9) == 6, 1, 'last')-100;
gamma = test_mat(index, 1);
alpha_star = test_mat(index, 2);
P_star = test_mat(index, 3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star-0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off
    
text(0, 1, '(d)', 'FontSize', 22)

%% Formation of Period-10 orbit

index = find(test_mat(:,9) == 10, 1, 'first');
gamma = test_mat(index, 1);
alpha_star = test_mat(index, 2);
P_star = test_mat(index, 3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star-0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off
    
text(0, 1, '(e)', 'FontSize', 22)

%% Period-10 orbit about to hit vertex

index = find(test_mat(:,9) == 10, 1, 'last');
gamma = test_mat(index, 1);
alpha_star = test_mat(index, 2);
P_star = test_mat(index, 3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star-0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off
    
text(0, 1, '(f)', 'FontSize', 22)

%% Period-14 orbit formation

index = find(test_mat(:,9) == 14, 1, 'first');
gamma = test_mat(index, 1);
alpha_star = test_mat(index, 2);
P_star = test_mat(index, 3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star-0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off
    
text(0, 1, '(g)', 'FontSize', 22)

%% Period-14 orbit about to hit vertex

index = find(test_mat(:,9) == 14, 1, 'last')-50;
gamma = test_mat(index, 1);
alpha_star = test_mat(index, 2);
P_star = test_mat(index, 3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star-0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off
    
text(0, 1, '(h)', 'FontSize', 22)




%% Points which support the period-N orbit in the Appendix

% For gamma = pi/2,

h = 1; gamma = pi/2;
P = linspace(0.001, 1+h/sin(gamma)-0.001, 10000); % all the possible initial positions on base and right adj side

row_index = 1;

gamma = test_mat(row_index,1);
alpha = test_mat(row_index,2);
period = test_mat(row_index,end);

for j=1:length(P)
    [side, angle, position] = parallelogram_map(h, gamma, alpha, P(j), period);
    angle_6(j) = angle(period+1); P_6(j) = position(period+1);
end

set(gca,'fontsize',140, 'FontWeight', 'Bold')

plot(P, P_6, 'ro', 'MarkerSize', 4); pbaspect([1 1 1]); hold on
plot(P, P, 'b-', 'LineWidth', 1)

xlabel('$P_0$','Interpreter','latex','String','{\boldmath$P_0$}'); 
ylabel('$P_2=g^{(2)}(\alpha_0, P_0)$', 'Interpreter','latex','String', '{\boldmath$P_2=g^{(2)}(\alpha_0, P_0)$}')
xticks([0 1 1+1/sin(gamma)]); xticklabels({'A', 'B', 'C'})
yticks([0 1 1+1/sin(gamma)]); yticklabels({'A', 'B', 'C'})

set(gca,'XLim',[0 1+1/sin(gamma)], 'YLim',[0 1+1/sin(gamma)], 'FontSize', 15);
%title(sprintf('gamma = %0.3f', gamma))
set(gcf,'color','w');

text(-0.4, 1+1/sin(gamma), '(a)', 'Color','k', 'FontSize', 22)


%% Period-2 about to hit vertex

row_index = find(test_mat(:,9) == 2, 1, 'last')-1000;

h = 1;
gamma = test_mat(row_index,1);
alpha = test_mat(row_index,2);
period = test_mat(row_index,end);

P = linspace(0.001, 1+h/sin(gamma)-0.001, 10000); % all the possible initial positions on base and right adj side

for j=1:length(P)
    [side, angle, position] = parallelogram_map(h, gamma, alpha, P(j), period);
    angle_6(j) = angle(period+1); P_6(j) = position(period+1);
end

set(gca,'fontsize',140, 'FontWeight', 'Bold')

plot(P, P, 'b-', 'LineWidth', 4); pbaspect([1 1 1]); hold on
plot(P, P_6, 'ro', 'MarkerSize', 4)

xlabel('$P_0$','Interpreter','latex','String','{\boldmath$P_0$}'); 
ylabel('$P_2=g^{(2)}(\alpha_0, P_0)$', 'Interpreter','latex','String', '{\boldmath$P_2=g^{(2)}(\alpha_0, P_0)$}')
xticks([0 1 1+1/sin(gamma)]); xticklabels({'A', 'B', 'C'})
yticks([0 1 1+1/sin(gamma)]); yticklabels({'A', 'B', 'C'})

set(gca,'XLim',[0 1+1/sin(gamma)], 'YLim',[0 1+1/sin(gamma)], 'FontSize', 15);
%title(sprintf('gamma = %0.3f', gamma))
set(gcf,'color','w');

text(-0.4, 1+1/sin(gamma), '(b)', 'Color','k', 'FontSize', 22)


%% Formation of Period-6 orbit

row_index = find(test_mat(:,9) == 6, 1, 'first');

h = 1;
gamma = test_mat(row_index,1);
alpha = test_mat(row_index,2);
period = test_mat(row_index,end);

P = linspace(0.001, 1+h/sin(gamma)-0.001, 10000); % all the possible initial positions on base and right adj side

for j=1:length(P)
    [side, angle, position] = parallelogram_map(h, gamma, alpha, P(j), period);
    angle_6(j) = angle(period+1); P_6(j) = position(period+1);
end

set(gca,'fontsize',140, 'FontWeight', 'Bold')

plot(P, P, 'b-', 'LineWidth', 4); pbaspect([1 1 1]); hold on
plot(P, P_6, 'ro', 'MarkerSize', 4)

xlabel('$P_0$','Interpreter','latex','String','{\boldmath$P_0$}'); 
ylabel('$P_6=g^{(6)}(\alpha_0, P_0)$', 'Interpreter','latex','String', '{\boldmath$P_6=g^{(6)}(\alpha_0, P_0)$}')
xticks([0 1 1+1/sin(gamma)]); xticklabels({'A', 'B', 'C'})
yticks([0 1 1+1/sin(gamma)]); yticklabels({'A', 'B', 'C'})

set(gca,'XLim',[0 1+1/sin(gamma)], 'YLim',[0 1+1/sin(gamma)], 'FontSize', 15);
%title(sprintf('gamma = %0.3f', gamma))
set(gcf,'color','w');

text(-0.4, 1+1/sin(gamma), '(c)', 'Color','k', 'FontSize', 22)


%% Period-6 orbit about to hit vertex

row_index = find(test_mat(:,9) == 6, 1, 'last')-100;

h = 1;
gamma = test_mat(row_index,1);
alpha = test_mat(row_index,2);
period = test_mat(row_index,end);

P = linspace(0.001, 1+h/sin(gamma)-0.001, 10000); % all the possible initial positions on base and right adj side

for j=1:length(P)
    [side, angle, position] = parallelogram_map(h, gamma, alpha, P(j), period);
    angle_6(j) = angle(period+1); P_6(j) = position(period+1);
end

set(gca,'fontsize',140, 'FontWeight', 'Bold')

plot(P, P, 'b-', 'LineWidth', 4); pbaspect([1 1 1]); hold on
plot(P, P_6, 'ro', 'MarkerSize', 4)

xlabel('$P_0$','Interpreter','latex','String','{\boldmath$P_0$}'); 
ylabel('$P_6=g^{(6)}(\alpha_0, P_0)$', 'Interpreter','latex','String', '{\boldmath$P_6=g^{(6)}(\alpha_0, P_0)$}')
xticks([0 1 1+1/sin(gamma)]); xticklabels({'A', 'B', 'C'})
yticks([0 1 1+1/sin(gamma)]); yticklabels({'A', 'B', 'C'})

set(gca,'XLim',[0 1+1/sin(gamma)], 'YLim',[0 1+1/sin(gamma)], 'FontSize', 15);
%title(sprintf('gamma = %0.3f', gamma))
set(gcf,'color','w');

text(-0.4, 1+1/sin(gamma), '(d)', 'Color','k', 'FontSize', 22)



%% Formation of Period-10 orbit

row_index = find(test_mat(:,9) == 10, 1, 'first');

h = 1;
gamma = test_mat(row_index,1);
alpha = test_mat(row_index,2);
period = test_mat(row_index,end);

P = linspace(0.001, 1+h/sin(gamma)-0.001, 10000); % all the possible initial positions on base and right adj side

for j=1:length(P)
    [side, angle, position] = parallelogram_map(h, gamma, alpha, P(j), period);
    angle_6(j) = angle(period+1); P_6(j) = position(period+1);
end

set(gca,'fontsize',140, 'FontWeight', 'Bold')

plot(P, P, 'b-', 'LineWidth', 4); pbaspect([1 1 1]); hold on
plot(P, P_6, 'ro', 'MarkerSize', 4)

xlabel('$P_0$','Interpreter','latex','String','{\boldmath$P_0$}'); 
ylabel('$P_{10}=g^{(10)}(\alpha_0, P_0)$', 'Interpreter','latex','String', '{\boldmath$P_{10}=g^{(10)}(\alpha_0, P_0)$}')
xticks([0 1 1+1/sin(gamma)]); xticklabels({'A', 'B', 'C'})
yticks([0 1 1+1/sin(gamma)]); yticklabels({'A', 'B', 'C'})

set(gca,'XLim',[0 1+1/sin(gamma)], 'YLim',[0 1+1/sin(gamma)], 'FontSize', 15);
%title(sprintf('gamma = %0.3f', gamma))
set(gcf,'color','w');

text(-0.4, 1+1/sin(gamma), '(e)', 'Color','k', 'FontSize', 22)


%% Period-10 orbit about to hit vertex

row_index = find(test_mat(:,9) == 10, 1, 'last');

h = 1;
gamma = test_mat(row_index,1);
alpha = test_mat(row_index,2);
period = test_mat(row_index,end);

P = linspace(0.001, 1+h/sin(gamma)-0.001, 10000); % all the possible initial positions on base and right adj side

for j=1:length(P)
    [side, angle, position] = parallelogram_map(h, gamma, alpha, P(j), period);
    angle_6(j) = angle(period+1); P_6(j) = position(period+1);
end

set(gca,'fontsize',140, 'FontWeight', 'Bold')

plot(P, P, 'b-', 'LineWidth', 4); pbaspect([1 1 1]); hold on
plot(P, P_6, 'ro', 'MarkerSize', 4)

xlabel('$P_0$','Interpreter','latex','String','{\boldmath$P_0$}'); 
ylabel('$P_{10}=g^{(10)}(\alpha_0, P_0)$', 'Interpreter','latex','String', '{\boldmath$P_{10}=g^{(10)}(\alpha_0, P_0)$}')
xticks([0 1 1+1/sin(gamma)]); xticklabels({'A', 'B', 'C'})
yticks([0 1 1+1/sin(gamma)]); yticklabels({'A', 'B', 'C'})

set(gca,'XLim',[0 1+1/sin(gamma)], 'YLim',[0 1+1/sin(gamma)], 'FontSize', 15);
%title(sprintf('gamma = %0.3f', gamma))
set(gcf,'color','w');

text(-0.4, 1+1/sin(gamma), '(f)', 'Color','k', 'FontSize', 22)


%% Formation of Period-14 orbit

row_index = find(test_mat(:,9) == 14, 1, 'first');

h = 1;
gamma = test_mat(row_index,1);
alpha = test_mat(row_index,2);
period = test_mat(row_index,end);

P = linspace(0.001, 1+h/sin(gamma)-0.001, 10000); % all the possible initial positions on base and right adj side

for j=1:length(P)
    [side, angle, position] = parallelogram_map(h, gamma, alpha, P(j), period);
    angle_6(j) = angle(period+1); P_6(j) = position(period+1);
end

set(gca,'fontsize',140, 'FontWeight', 'Bold')

plot(P, P, 'b-', 'LineWidth', 4); pbaspect([1 1 1]); hold on
plot(P, P_6, 'ro', 'MarkerSize', 4)

xlabel('$P_0$','Interpreter','latex','String','{\boldmath$P_0$}'); 
ylabel('$P_{14}=g^{(14)}(\alpha_0, P_0)$', 'Interpreter','latex','String', '{\boldmath$P_{14}=g^{(14)}(\alpha_0, P_0)$}')
xticks([0 1 1+1/sin(gamma)]); xticklabels({'A', 'B', 'C'})
yticks([0 1 1+1/sin(gamma)]); yticklabels({'A', 'B', 'C'})

set(gca,'XLim',[0 1+1/sin(gamma)], 'YLim',[0 1+1/sin(gamma)], 'FontSize', 15);
%title(sprintf('gamma = %0.3f', gamma))
set(gcf,'color','w');

text(-0.4, 1+1/sin(gamma), '(g)', 'Color','k', 'FontSize', 22)


%% Period-14 orbit about to hit vertex

row_index = find(test_mat(:,9) == 14, 1, 'last')-50;

h = 1;
gamma = test_mat(row_index,1);
alpha = test_mat(row_index,2);
period = test_mat(row_index,end);

P = linspace(0.001, 1+h/sin(gamma)-0.001, 10000); % all the possible initial positions on base and right adj side

for j=1:length(P)
    [side, angle, position] = parallelogram_map(h, gamma, alpha, P(j), period);
    angle_6(j) = angle(period+1); P_6(j) = position(period+1);
end

set(gca,'fontsize',140, 'FontWeight', 'Bold')

plot(P, P, 'b-', 'LineWidth', 4); pbaspect([1 1 1]); hold on
plot(P, P_6, 'ro', 'MarkerSize', 4)

xlabel('$P_0$','Interpreter','latex','String','{\boldmath$P_0$}'); 
ylabel('$P_{14}=g^{(14)}(\alpha_0, P_0)$', 'Interpreter','latex','String', '{\boldmath$P_{14}=g^{(14)}(\alpha_0, P_0)$}')
xticks([0 1 1+1/sin(gamma)]); xticklabels({'A', 'B', 'C'})
yticks([0 1 1+1/sin(gamma)]); yticklabels({'A', 'B', 'C'})

set(gca,'XLim',[0 1+1/sin(gamma)], 'YLim',[0 1+1/sin(gamma)], 'FontSize', 15);
%title(sprintf('gamma = %0.3f', gamma))
set(gcf,'color','w');

text(-0.4, 1+1/sin(gamma), '(h)', 'Color','k', 'FontSize', 22)






%% All the different angles and positions visited 

Nmax = test_mat(end,9);
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
zlabel('P - position', 'FontSize', 22)

xticks([0 pi/8 pi/4 3*pi/8 pi/2])
xticklabels({'0', '\pi/8', '\pi/4', '3\pi/8', '\pi/2'})

period = prod(6:4:Nmax);
title(sprintf('Continuation Period-%d', period))

[x y] = meshgrid(0:0.01:pi/2, 0:0.01:pi); 
z0 = ones(size(x));
z1 = 1 + h./sin(x);    
z2 = 2 + h./sin(x);
z3 = 2 + 2*h./sin(x);

% For dissertation
x = linspace(0, pi/2, 100);
plot3(x, zeros(size(x)), zeros(size(x)), 'r-', 'LineWidth', 2); hold on
plot3(x, zeros(size(x)), ones(size(x)), 'r-', 'LineWidth', 2); hold on
plot3(x, zeros(size(x)), 1 + h./sin(x), 'r-', 'LineWidth', 2); hold on
plot3(x, zeros(size(x)), 2 + h./sin(x), 'r-', 'LineWidth', 2); hold on
plot3(x, zeros(size(x)), 2 + 2*h./sin(x), 'r-', 'LineWidth', 2); hold on

view(0, 0) % az, el

text(-0.1, 0, 8, '(a)', 'FontSize', 22)
set(gcf,'color','w'); % background color to white

%% Plot of gamma against the period

plot(test_mat(:,1), test_mat(:,9), 'o')
xlabel('\gamma'); ylabel('Period')
xticks([0 pi/8 pi/4 3*pi/8 pi/2])
xticklabels({'0', '\pi/8', '\pi/4', '3\pi/8', '\pi/2'})

set(gca,'XLim',[0 pi/2]);
set(gca,'FontSize', 20)
set(gcf,'color','w'); % background color to white

text(-0.15, 60, '(b)', 'FontSize', 22)

%% Newton Steps (in Appendix)

% Column 1 is gamma
% Column 2 is initial angle, alpha
% Column 3 is initial position, P
% Column 4 is the area
% Column 5 is the number of Newton steps
% Column 6 & 7 are function evaluations for g^(N)(alpha,P)-P and dA/dP
% Column 8 is the function evaluation for alpha, f^(N)(alpha,P)-alpha

plot(test_mat(:,1),test_mat(:,5), 'o', 'Markersize', 4)
xlabel('$\gamma$', 'FontSize', 22); ylabel('Number of Newton Steps', 'FontSize', 22)

set(gca,'XLim',[0 pi/2], 'YLim', [0 max(test_mat(:,5))], 'FontSize', 15)
xticks([0 pi/8 pi/4 3*pi/8 pi/2])
xticklabels({'0', '\pi/8', '\pi/4', '3\pi/8', '\pi/2'})
yticks([0 1 2 3])
yticklabels({'0', '1', '2', '3'})
grid on
set(gcf,'color','w'); % background color to white
text(-0.15, 3, '(a)', 'FontSize', 22)

%% Period-adding area (in Appendix)

area = test_mat(:,4);

plot(test_mat(:,1),area, 'o', 'Markersize', 4)
xlabel('$\gamma$', 'FontSize', 22); ylabel('Area', 'FontSize', 22)

set(gca,'XLim',[0 pi/2], 'YLim', [min(area)-0.1 max(area)], 'FontSize', 15)
xticks([0 pi/8 pi/4 3*pi/8 pi/2])
xticklabels({'0', '\pi/8', '\pi/4', '3\pi/8', '\pi/2'})
yticks([0.9 1 1.1 1.2])
yticklabels({'0.9', '1', '1.1', '1.2'})
grid on
set(gcf,'color','w'); % background color to white
text(-0.15, 1.2, '(b)', 'FontSize', 22)



%% Introduction to our square map.
% Fig 3.1
N = 4; % number of iterations

gamma = pi/2; h = 1;   % A square
alpha_star = pi/4-0.4;
P_star = 0.2;

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold on

% Need to relabel the corners as 0,1,2,3

% Plotting dots for each bounce, labels and angle arcs
for i = 1:N
    if side(i) == 0
        plot(position(i), 0, 'go', 'MarkerFaceColor', 'g')
        text(position(i), 0-0.05, 'P_{0}', 'FontSize', 18, 'Color', 'g')
        ang([position(i) 0], 0.2, [0 alpha_star],'r-'); % angle
        text(position(i)+0.25, 0.05, '\alpha_{0}', 'Color', 'r', 'FontSize', 18)
    elseif side(i) == 1
        plot(1, position(i) - 1, 'go', 'MarkerFaceColor', 'g')
        text(1+0.03, position(i) - 1, 'P_{1}', 'FontSize', 18, 'Color', 'g')
        ang([1 position(i)-1], 0.1, [0 pi-alpha_star],'r-'); % angle
        text(1-0.1, position(i)-1+0.15, '\alpha_{1}', 'Color', 'r', 'FontSize', 18)
    elseif side(i) == 2
        plot(-position(i) + 3, 1, 'go', 'MarkerFaceColor', 'g')
        text(-position(i) + 3, 1+0.05, 'P_{3}', 'FontSize', 18, 'Color', 'g')
        ang([-position(i) + 3 1], 0.1, [pi 2*pi-alpha_star],'r-'); % angle
        text(-position(i) + 3, 1-0.15, '\alpha_{3}', 'Color', 'r', 'FontSize', 18)
    elseif side(i) == 3
        plot(0, -position(i) + 4, 'go', 'MarkerFaceColor', 'g')
        text(0-0.08, -position(i) + 4, 'P_{2}', 'FontSize', 18, 'Color', 'g')
        ang([0 -position(i) + 4], 0.1, [pi 2*pi+alpha_star],'r-'); % angle
        text(0+0.1, -position(i) + 4-0.1, '\alpha_{2}', 'Color', 'r', 'FontSize', 18)
    end
end
    
%plot()


%% Bifucation diagram for period-6 orbit termination

h = 1;
gamma = pi/2;
period = 6;

init_alpha_1 = atan(2); init_P_1 = 0.25;
init_alpha_2 = atan(1/2); init_P_2 = 0.5;

bigmat_0p25 = parallelogram_newton_solver(h, gamma, init_alpha_1, init_P_1, period);
bigmat_0p5 = parallelogram_newton_solver(h, gamma, init_alpha_2, init_P_2, period);

for j = 1:10:(size(bigmat_0p25, 1)-1)    % last entry is the one we terminate at, so is scuffed
    gamma = bigmat_0p25(j, 1);
    alpha = bigmat_0p25(j, 2); P = bigmat_0p25(j, 3);
    
    [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
    
    plot3(gamma * ones(size(a)), a, p, 'bo', 'MarkerSize', 1); hold on
    %plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'bo', 'MarkerSize', 1); hold on
end
view(0,0)


for j = 1:10:(size(bigmat_0p5, 1)-1)    % last entry is the one we terminate at, so is scuffed
    gamma = bigmat_0p5(j, 1);
    alpha = bigmat_0p5(j, 2); P = bigmat_0p5(j, 3);
    
    [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
    
    plot3(gamma * ones(size(a)), a, p, 'go', 'MarkerSize', 1); hold on
    %plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'bo', 'MarkerSize', 1); hold on
end
view(0,0)

x = linspace(1.35, pi/2, 100);
plot3(x, zeros(size(x)), zeros(size(x)), 'r-', 'LineWidth', 2); hold on
plot3(x, zeros(size(x)), ones(size(x)), 'r-', 'LineWidth', 2); hold on
plot3(x, zeros(size(x)), 1 + h./sin(x), 'r-', 'LineWidth', 2); hold on
plot3(x, zeros(size(x)), 2 + h./sin(x), 'r-', 'LineWidth', 2); hold on
plot3(x, zeros(size(x)), 2 + 2*h./sin(x), 'r-', 'LineWidth', 2); hold on

xlabel('\gamma'); zlabel('Positions visited')
xticks([1.35 1.46 pi/2])
xticklabels({'1.35', '1.46', '\pi/2'})
zticks([0 1 2 3 4])
zticklabels({'0', '1', '2', '3' ,'4'})
set(gcf,'color','w'); % background color to white
set(gca,'FontSize', 20,'XLim',[1.35 pi/2]); box on

legend({'\alpha_0 = 0.25, P_0 = atan(2)', '\alpha_0=0.5, P_0 = atan(1/2)', 'Vertex position'}, 'FontSize',14)
%legend boxoff 


%% Parallelogram used in cutting sequence geometric proof Chapter 4

h = 0.7; % Same parallelogram
gamma = pi/4+0.2;

% First trajectory
alpha_star = 0.3;
P_star = 0.5;
N = 1;
[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold on
%plot(P_star, 0, 'go', 'MarkerFaceColor', 'g') % A nice dot
text(P_star, -0.05, '$P_{0}$', 'Color', 'g', 'FontSize', 18)
ang([P_star 0], 0.15, [0 alpha_star],'b-'); % angle
ang([P_star 0], 0.1, [pi pi-2],'b-'); % angle
text(P_star+0.17, 0.03, '$\alpha$', 'Color', 'b', 'FontSize', 18)

alpha_star = 1.9; P_star = position(2);
[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold on
text(1.15, 0.19, '$P_{1}$', 'Color', 'g', 'FontSize', 18)
%ang([1.15 0.19], 0.08, [1.2 pi-0.24],'b-'); % angle
ang([1.15 0.19], 0.08, [0 2*pi],'b-'); % angle
text(1.15-0.05, 0.3, '$\beta$', 'Color', 'b', 'FontSize', 18)

alpha_star = pi-0.5; P_star = position(2);
[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold on
text(0.22, 0.42, '$P_{2}$', 'Color', 'g', 'FontSize', 18)
ang([0.28 0.42], 0.08, [0 2*pi],'b-'); % angle
text(0.3, 0.32, '$\epsilon$', 'Color', 'b', 'FontSize', 18)

alpha_star = 1.15; P_star = position(2);
[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold on
text(0.8, 0.73, '$P_{3}$', 'Color', 'g', 'FontSize', 18)
ang([0.8 0.73], 0.08, [0 2*pi],'b-'); % angle
text(0.88, 0.65, '$\theta$', 'Color', 'b', 'FontSize', 18)

% Angle arcs
% Bottom left
ang([0 0], 0.1, [0 gamma],'r-'); % angle
text(0.1, 0.07, '$\gamma$', 'Color', 'r', 'FontSize', 18)
% Bottom right
ang([1 0], 0.05, [pi/2-0.57 pi],'r-'); % angle
text(0.9, 0.07, '$\pi - \gamma$', 'Color', 'r', 'FontSize', 18)
% Top left
ang([h/tan(gamma) h], 0.05, [pi+gamma 2*pi],'r-'); % angle
text(h/tan(gamma)+0.03, h-0.05, '$\pi - \gamma$', 'Color', 'r', 'FontSize', 18)
% Top right
ang([1+h/tan(gamma) h], 0.1, [pi+gamma pi],'r-'); % angle
text(1+h/tan(gamma)-0.17, h-0.05, '$\pi - \gamma$', 'Color', 'r', 'FontSize', 18)

% Middle dot
plot(0.643, 0.323, 'ro', 'MarkerFaceColor', 'r')
text(0.68, 0.35, 'O', 'Color', 'r', 'FontSize', 18)


%% alpha_m and alpha_m' values in the Appendix

x = linspace(-pi, pi, 100);    % alpha_m
gamma = 0.6;

numerator = sin(gamma + x);
denominator = sin(x);

plot(x, numerator, 'b','LineWidth',2); hold on
plot(x, denominator, 'r','LineWidth',2); hold on
plot(x, zeros(size(x)), 'k') % y axis
set(gca,'XLim',[0 pi],'YLim',[-1 1], 'FontSize', 15)
set(gcf,'color','w'); % background color to white

xticks([0 pi-gamma pi]);
xticklabels({'0', '$\pi-\gamma$', '$\pi$'}); xlabel('$\alpha_m$', 'FontSize', 15)
yticks([-1 0 1]);
yticklabels({'-1', '0', '1'});

legend({'sin$\left(\gamma+\alpha_m\right)$', 'sin$\left(\alpha_m\right)$'}, 'FontSize',14, 'Interpreter','latex')
legend boxoff   


%% Period-2 example trajectory in the Appendix

h = 1; gamma = pi/2-0.4;

alpha_star = pi/2; P_star = 1.4; 
N = 10;

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold on

% starting just off vertex 1
[side alpha position] = parallelogram_map(h, gamma, alpha_star, 1.01, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold on

% starting just off vertex 3
[side alpha position] = parallelogram_map(h, gamma, alpha_star, 2+h/sin(gamma)+0.01, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold on
