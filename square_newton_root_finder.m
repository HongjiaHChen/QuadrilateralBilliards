
init_alpha = pi/3; init_P = 0.1;
init_guess = [init_alpha; init_P]; % init_angle and init_pos

period = 4; % what period are we seaching for?


%[F_alpha, F_P] = square_map(init_alpha, init_P, period);
%F_n = [F_alpha(period+1); F_P(period+1)];

%next_guess = init_guess - inv(square_jacobian(init_alpha, init_P, period) - eye(2)) * (F_n - init_guess);

alpha = init_alpha; P = init_P;
guess=  init_guess;

for j=1:10
    % F^{N}(\alpha, \P) term
   [F_alpha, F_P] = square_map(alpha, P, period);
    F_n = [F_alpha(period+1); F_P(period+1)];
    
   guess = guess - inv(square_jacobian(alpha, P, period) - eye(2)) * (F_n - guess);
   
   alpha = guess(1); P = guess(2);
end