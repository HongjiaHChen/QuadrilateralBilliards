function bigmat = parallelogram_newton_solver(h, gamma, init_alpha, init_P, period)

% This function keeps h constant and varies gamma as we track a periodic
% orbit for varying gamma.
% User arguments
%gamma: angle in parallelogram
%h: height of parallelogram
%init_alpha: inital angle which gives a periodic orbit with maximal area in
%the parallelogram with angle gamma and height h
%init_P: initial position as above
%period:  what period are we seaching for?

eps = 1e-6;   % for finite differencing

gamma_perturb = 1e-5; %1e-3; % perturbing gamma of the parallelogram

gammas = [gamma];             % The gammas, filled in for square already
alpha_P_mat = [init_alpha init_P]; % This alpha and P is the sol. for the square. 

areas = square_area_V2(alpha_P_mat(1), alpha_P_mat(2), period);   % store the areas as well
newton_steps = [0];                        % number of steps to convergence
function_eval = zeros(1,2);                       % what the function evaluates to (should be close to 0).

for i=1:25000   % gradually change the gamma of the parallelogram
    
    gammas(i+1) = gammas(1)-gamma_perturb*i;

    gamma = gammas(i+1); % gamma of parallelogram
    
    % The initial guess for alpha and P depends on the last non-NaN values
    % From https://au.mathworks.com/matlabcentral/answers/29481-last-non-nan-observation
     B = ~isnan(alpha_P_mat);
     % indices
     Indices = arrayfun(@(x) find(B(:, x), 1, 'last'), 1:size(alpha_P_mat, 2));
     % values
     Values = arrayfun(@(x,y) alpha_P_mat(x,y), Indices, 1:size(alpha_P_mat, 2));
     
    guesses = Values; % First column: init_angle and init_pos  
    
    try % the parallelogram_jacobian_area func causes some errors as no matter what step we take, 
        % we seem to run into a vertex. this makes calculations of the
        % derivatives using finite differencing difficult
        
        max_newton_iter = 4;
        for j=1:max_newton_iter
            guess = guesses(j,:)';  % column vector
            alpha = guess(1); P = guess(2);

            % F^{N}(\alpha, \P) term. Where are we after a period number of bounces
            [F_side F_alpha, F_P] = parallelogram_map(h, gamma, alpha, P, period);
            F_n = [F_alpha(period+1); F_P(period+1)];

            dA_dP = (parallelogram_area(h, gamma, alpha, P+eps, period)-parallelogram_area(h, gamma, alpha, P-eps, period))/(2*eps);

            % if our initial guess happens satisfy the root finding problem
            % as well as the angle returning to be the same
            % IMPT: will always be the case for period-2 orbits as area
            % remains constant
            
            
            
            %[abs(F_n(2)-P); abs(dA_dP); abs(F_n(1)-alpha)]'
            
            
            
            if (j == 1) && (abs(F_n(2)-P) < 1e-6) && (abs(dA_dP) < 1e-5) && (abs(F_n(1)-alpha) < 1e-6)
                areas(i+1) = parallelogram_area(h, gamma, alpha, P, period);
                newton_steps(i+1) = 0;    % our initial guess was correct
                function_eval(i+1,:) = [F_n(2)-P dA_dP];
                alpha_eval(i+1) = F_n(1)-alpha;
                break
            end
            
            % otherwise do Newton's method
            guess = guess - inv(parallelogram_jacobian_area(h, gamma, alpha, P, period)) * [F_n(2)-P; dA_dP];

            alpha = guess(1); P = guess(2);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TEMP SOLUTION for when we hop out of range
            %alpha = mod(alpha, pi); P = mod(P, 1);
            alpha = mod(alpha, pi); P = mod(P, 2+2*h/sin(gamma));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            guesses(j+1,:) = [alpha P];
            
            if (norm(guesses(j+1,:)- guesses(j,:))/norm(guesses(j,:)) < 1e-6)% relative error
                
                % For function evaluation, have we reached [0,0] yet?
                [F_side F_alpha, F_P] = parallelogram_map(h, gamma, alpha, P, period);
                F_n = [F_alpha(period+1); F_P(period+1)];
                dA_dP = (parallelogram_area(h, gamma, alpha, P+eps, period)-parallelogram_area(h, gamma, alpha, P-eps, period))/(2*eps);
                
                %[(abs(F_n(2)-guesses(j+1,2))) (abs(dA_dP))
                %abs(F_n(1)-guesses(j+1,1))] % debugging
                
                if (abs(F_n(2)-guesses(j+1,2)) < 1e-6) && (abs(dA_dP)< 1e-8) && (abs(F_n(1)-guesses(j+1,1)) < 1e-6) %angle must be the same
                    areas(i+1) = parallelogram_area(h, gamma, alpha, P, period);
                    newton_steps(i+1) = j;    % we converged on jth iteration
                    function_eval(i+1,:) = [F_n(2)-guesses(j+1,2) dA_dP];
                    alpha_eval(i+1) = F_n(1)-guesses(j+1,1);
                    break
                end
                
            end
        end
        
    if j == max_newton_iter   %we reached end of Newton's and didn't converge
        
        %guesses(j+1,:) = [NaN NaN];
        alpha_P_mat(i+1,:) = guesses(j+1,:);
        areas(i+1) = parallelogram_area(h, gamma, alpha, P, period);     %NaN;
        newton_steps(i+1) = j;                                           %NaN;    
        function_eval(i+1,:) = F_n;                                      %[NaN NaN];

        alpha_eval(i+1) = F_n(1)-guesses(j+1,1);                         %NaN;
        
        break % terminate for loop as something has gone wrong (eg: continuation makes trajectory hits vertex)
    end
    
    
    alpha_P_mat(i+1,:) = guesses(end,:);  % the converged value
    
    
    catch
        % Did not converge to anything as an error occured
        guesses(j+1,:) = [NaN NaN];
        alpha_P_mat(i+1,:) = guesses(j+1,:);
        areas(i+1) = NaN;
        newton_steps(i+1) = NaN;    
        function_eval(i+1,:) = [NaN NaN];

        alpha_eval(i+1) = NaN;
        
        break
        %continue; % Jump to next iteration
        
    end
    
    
end

bigmat = [gammas' alpha_P_mat areas' newton_steps' function_eval alpha_eval'];


end

