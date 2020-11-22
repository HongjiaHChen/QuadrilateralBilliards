function J = parallelogram_jacobian_area(h, gamma, alpha_star, P_star, period)
% Function computes the Jacobian at a certain point using central
% differencing for the derivatives.

eps = 1e-6; % For the finite difference. We use the central finite difference method

% How many iteations should we go ahead? (How many compositions of the
% maps?)
N = period;

try % central difference approx.
    % dg/d(alpha) = (g(alpha + eps, P) - g(alpha - eps, P)) / (2*eps)      
    % g is equation for position
    [s dg_da_big_alpha dg_da_big_P] = parallelogram_map(h, gamma, alpha_star+eps, P_star, period);       % g(alpha+eps, P)
    [s dg_da_small_alpha dg_da_small_P] = parallelogram_map(h, gamma, alpha_star-eps, P_star, period);   % g(alpha-eps, P)
    dg_da = (dg_da_big_P(period+1)-dg_da_small_P(period+1))/(2*eps);               % approximately dg/d(alpha)
    
catch % if central difference fails (i.e: we hit a vertex), we use forward difference approx.
    % dg/d(alpha) = (g(alpha + eps, P) - g(alpha, P)) / (eps)   
    
    try % forward approx
    [s dg_da_big_alpha dg_da_big_P] = parallelogram_map(h, gamma, alpha_star+eps, P_star, period);       % g(alpha+eps, P)
    [s dg_da_small_alpha dg_da_small_P] = parallelogram_map(h, gamma, alpha_star, P_star, period);   % g(alpha, P)
    dg_da = (dg_da_big_P(period+1)-dg_da_small_P(period+1))/(eps);               % approximately dg/d(alpha)
    
    
    catch % backward approx
        [s dg_da_big_alpha dg_da_big_P] = parallelogram_map(h, gamma, alpha_star, P_star, period);       % g(alpha, P)
        [s dg_da_small_alpha dg_da_small_P] = parallelogram_map(h, gamma, alpha_star-eps, P_star, period);   % g(alpha-eps, P)
        dg_da = (dg_da_big_P(period+1)-dg_da_small_P(period+1))/(eps);               % approximately dg/d(alpha)
    end
    
end


try % central difference approx.
    % dg/dP = (g(alpha, P + eps) - g(alpha, P - eps)) / (2*eps)
    [s dg_dP_big_alpha dg_dP_big_P] = parallelogram_map(h, gamma,alpha_star, P_star+eps, period);       % g(alpha, P+eps)
    [s dg_dP_small_alpha dg_dP_small_P] = parallelogram_map(h, gamma,alpha_star, P_star-eps, period);   % g(alpha, P-eps)
    dg_dP = (dg_dP_big_P(period+1)-dg_dP_small_P(period+1))/(2*eps);          % approximately dg/dP
    
catch % If central diff fails, we try forward approx.
    
    try % foward approx.
        % dg/dP = (g(alpha, P + eps) - g(alpha, P)) / (eps)
        [s dg_dP_big_alpha dg_dP_big_P] = parallelogram_map(h, gamma,alpha_star, P_star+eps, period);       % g(alpha, P+eps)
        [s dg_dP_small_alpha dg_dP_small_P] = parallelogram_map(h, gamma,alpha_star, P_star, period);   % g(alpha, P)
        dg_dP = (dg_dP_big_P(period+1)-dg_dP_small_P(period+1))/(eps);          % approximately dg/dP
        
    catch % if foward fails, we try backward approx
        [s dg_dP_big_alpha dg_dP_big_P] = parallelogram_map(h, gamma,alpha_star, P_star, period);       % g(alpha, P)
        [s dg_dP_small_alpha dg_dP_small_P] = parallelogram_map(h, gamma,alpha_star, P_star-eps, period);   % g(alpha, P-eps)
        dg_dP = (dg_dP_big_P(period+1)-dg_dP_small_P(period+1))/(eps);          % approximately dg/dP
        
    end
    
end

% d^2A/d(alpha)dP
try % Central difference
    d2A_dadP = (parallelogram_area(h, gamma,alpha_star+eps, P_star+eps, period)-parallelogram_area(h, gamma,alpha_star+eps, P_star-eps, period)-...
           parallelogram_area(h, gamma,alpha_star-eps, P_star+eps, period)+parallelogram_area(h, gamma,alpha_star-eps, P_star-eps, period))/(4*eps^2);

catch % Forward difference
    
    try  % Forward difference
        d2A_dadP = (parallelogram_area(h, gamma,alpha_star+eps, P_star+eps, period)-parallelogram_area(h, gamma,alpha_star+eps, P_star, period)-...
               parallelogram_area(h, gamma,alpha_star, P_star+eps, period)+parallelogram_area(h, gamma,alpha_star, P_star, period))/(4*eps^2);
    
    catch % Backward difference
        d2A_dadP = (parallelogram_area(h, gamma,alpha_star, P_star, period)-parallelogram_area(h, gamma,alpha_star, P_star-eps, period)-...
           parallelogram_area(h, gamma,alpha_star-eps, P_star, period)+parallelogram_area(h, gamma,alpha_star-eps, P_star-eps, period))/(4*eps^2);
    end
end


% d^2A/dP^2

try % Central difference
    d2A_dP2 = (parallelogram_area(h, gamma,alpha_star, P_star+eps, period)-2*parallelogram_area(h, gamma,alpha_star, P_star, period)+...
            parallelogram_area(h, gamma,alpha_star, P_star-eps, period))/(eps^2);
catch % If central fails
    try % Forward difference
        
        d2A_dP2 = (parallelogram_area(h, gamma,alpha_star, P_star+2*eps, period)-2*parallelogram_area(h, gamma,alpha_star, P_star+eps, period)+...
                parallelogram_area(h, gamma,alpha_star, P_star, period))/(eps^2);
            
    catch % Backward difference
        
        d2A_dP2 = (parallelogram_area(h, gamma,alpha_star, P_star, period)-2*parallelogram_area(h, gamma,alpha_star, P_star-eps, period)+...
                parallelogram_area(h, gamma,alpha_star, P_star-2*eps, period))/(eps^2);
    end
end


J = [dg_da dg_dP-1; d2A_dadP d2A_dP2];
% eig(J)

end

