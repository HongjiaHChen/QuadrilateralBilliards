function J = rec_jacobian_area(L_H, L_V, alpha_star, P_star, period)
% Function computes the Jacobian at a certain point using central
% differencing for the derivatives.

eps = 1e-6; % For the finite difference. We use the central finite difference method

% How many iteations should we go ahead? (How many compositions of the
% maps?)
N = period;


% dg/d(alpha) = (g(alpha + eps, P) - g(alpha - eps, P)) / (2*eps)      
% g is equation for position
[s dg_da_big_alpha dg_da_big_P] = rec_map(L_H, L_V, alpha_star+eps, P_star, period);       % g(alpha+eps, P)
[s dg_da_small_alpha dg_da_small_P] = rec_map(L_H, L_V, alpha_star-eps, P_star, period);   % g(alpha-eps, P)
dg_da = (dg_da_big_P(period+1)-dg_da_small_P(period+1))/(2*eps);               % approximately dg/d(alpha)

% dg/dP = (g(alpha, P + eps) - g(alpha, P - eps)) / (2*eps)
[s dg_dP_big_alpha dg_dP_big_P] = rec_map(L_H, L_V,alpha_star, P_star+eps, period);       % g(alpha, P+eps)
[s dg_dP_small_alpha dg_dP_small_P] = rec_map(L_H, L_V,alpha_star, P_star-eps, period);   % g(alpha, P-eps)
dg_dP = (dg_dP_big_P(period+1)-dg_dP_small_P(period+1))/(2*eps);          % approximately dg/dP

% d^2A/d(alpha)dP
d2A_dadP = (rec_area(L_H, L_V,alpha_star+eps, P_star+eps, period)-rec_area(L_H, L_V,alpha_star+eps, P_star-eps, period)-...
           rec_area(L_H, L_V,alpha_star-eps, P_star+eps, period)+rec_area(L_H, L_V,alpha_star-eps, P_star-eps, period))/(4*eps^2);

% d^2A/dP^2
d2A_dP2 = (rec_area(L_H, L_V,alpha_star, P_star+eps, period)-2*rec_area(L_H, L_V,alpha_star, P_star, period)+...
            rec_area(L_H, L_V,alpha_star, P_star-eps, period))/(eps^2);
        
J = [dg_da dg_dP-1; d2A_dadP d2A_dP2];
% eig(J)

end

