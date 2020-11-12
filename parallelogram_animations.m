% Note we need to restrict gamma to be in the interval (0, pi/2].

%% Testing

h = 3.487;  % height of parallelogram
gamma = pi/2; % bottom left angle of parallelogram


N = 10;  % number of iterations
position0 = 0.3659; alpha0 = atan(2); % user input

[s, a, p] = parallelogram_map(h, gamma, alpha0, position0, N)



%% Period-6 orbit in square
h = 1; gamma = pi/2; % Square
N = 10;
alpha0 = atan(2); position0 = 0.2; 

[side alpha position] = parallelogram_map(h, gamma, alpha0, position0, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position)


%% Rectangle
h = 1.5; gamma = pi/2; 
N = 10;
alpha0 = atan(2); position0 = 0.2; 

[side alpha position] = parallelogram_map(h, gamma, alpha0, position0, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position)

%% Random Parallelogram

h = 1.5; gamma = pi/2-0.1; 
N = 10;
alpha0 = atan(2); position0 = 0.2; 

[side alpha position] = parallelogram_map(h, gamma, alpha0, position0, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position)


%% Debug testing for parallelogram

h = 1; gamma = 0.4; 
N = 10;
alpha0 = pi/2; position0 = 1.8; 

[side alpha position] = parallelogram_map(h, gamma, alpha0, position0, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);


%% Very angled parallelogram testing

h = 1; gamma = 0.4; 
N = 10;
alpha0 = pi/2; position0 = 0.05; 

[side alpha position] = parallelogram_map(h, gamma, alpha0, position0, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);








%%

h = 3.5; gamma = pi/2-0.01;

alpha0 = atan(h*3/2); position0 = 0;

N = 10;

[side alpha position] = parallelogram_map(h, gamma, alpha0, position0, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);





%%

row = 4;

N = 50000;
h = 1; gamma = bigmat(row, 1);

alpha_star = bigmat(row, 2);
P_star = bigmat(row, 3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);





%% Using generalized diagonals. Intentionally shooting our billiard into the vertex.

N = 100;
h = 1; gamma = pi/2-0.2;

%alpha_star = 0.149; P_star = 0;
alpha_star =atan(h/(1+h/tan(gamma))); P_star = 0.1;   % Period 8 orbit!
% Doesn't work for other h though :(

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);



%% Using generalized diagonals. Intentionally shooting our billiard into the vertex.

N = 100;
h = 1; gamma = pi/2-0.2;

alpha_star =mod(-gamma+atan(-h/(1-h/tan(gamma))),pi); P_star = 0.5;   % this one goes from bottom right to top left. 
% doesn't quite work out

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);

%% This doesn't work as well
N = 10;
h = 1; gamma = pi/2-0.2;

alpha_star =gamma; P_star = 0.5;   

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);



%% Looping through to produce images like Franks

N = 30; % how long the billiard goes for on each iteration
N_iter = 5000; % number of iterations

h = 1; % height of parallelogram

gammas_vec = linspace(0.3, pi/2, N_iter); % gammas for each iteration

init_pos = 0.5;

for j=1:N_iter  % change the heights with each iteration
    gamma = gammas_vec(j); 
    init_angle = atan(h/(1+h/tan(gamma)));
    
    [side alpha position] = parallelogram_map(h, gamma, init_angle, init_pos, N);
    P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off
    
    Image = getframe(gcf);
    folder = 'C:/Users/Work/Desktop/UoA/Year 4/Honours Project/Mathematical Billiards/MATLAB Code/Images/Parallelogram/Animation/30 iterations init_pos=0.5 Period 4';
    baseFileName = sprintf('test%d.jpg', j);
    fullFileName = fullfile(folder, baseFileName);
    imwrite(Image.cdata, fullFileName);
end



%%















%%
h = 1;
gamma = pi/2;
period = 6*10;

init_alpha_1 = atan(2/1);

bigmat_0p25 = parallelogram_newton_solver(h, gamma, init_alpha_1, 0.25, period)


%%
row = size(bigmat_0p25,1) - 1;
bigmat_0p25(row,:);

N = 8;

h = 1; gamma = bigmat_0p25(row, 1) ;

alpha_star = bigmat_0p25(row, 2); P_star = bigmat_0p25(row, 3);   

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);

%%
for j=1:size(bigmat_0p25,1)  % iterating over each row in the matrix
    
    if ~isnan(bigmat_0p25(j, 2))
        N = 50;

        h = 1; gamma = bigmat_0p25(j, 1) ; %1.5608;

        alpha_star = bigmat_0p25(j, 2); P_star = bigmat_0p25(j, 3);   

        [side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
        P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

        Image = getframe(gcf);
%         folder = 'C:/Users/Work/Desktop/UoA/Year 4/Honours Project/Mathematical Billiards/MATLAB Code/Images/Parallelogram/Animation/50 iterations follow Period 6';
%         baseFileName = sprintf('test%d.jpg', j);
%         fullFileName = fullfile(folder, baseFileName);
%         imwrite(Image.cdata, fullFileName);
    end
    
end













%%
row = size(bigmat,1)-1;
bigmat(row,:)

%% Testing, feel free to change this
N = 8;

h = 1; gamma = bigmat(row, 1) ; %1.5608;

alpha_star = bigmat(row, 2); P_star = bigmat(row, 3)+0.01;   

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);



%%

for j=1:size(bigmat,1)  % iterating over each row in the matrix
    
    if ~isnan(bigmat(j, 2))
        N = 50;

        h = 1; gamma = bigmat(j, 1) ; %1.5608;

        alpha_star = bigmat(j, 2); P_star = bigmat(j, 3);   

        [side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
        P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

        Image = getframe(gcf);
        folder = 'C:/Users/Work/Desktop/UoA/Year 4/Honours Project/Mathematical Billiards/MATLAB Code/Images/Parallelogram/Animation/50 iterations follow Period 6';
        baseFileName = sprintf('test%d.jpg', j);
        fullFileName = fullfile(folder, baseFileName);
        imwrite(Image.cdata, fullFileName);
    end
    
end







%% Parallelogram for dissertation
N = 0;
h = 1; gamma = pi/2-0.4;

alpha_star =gamma; P_star = 0.5;   

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);


% Angle Arcs
ang([0 0], 0.15, [0 gamma],'r-')
text(0.17, 0.17, '\gamma', 'Color', 'r', 'FontSize', 22)

text(1.3, 0.5, '$h$', 'Color', 'r', 'FontSize', 22,'Interpreter','latex')




