function J = parallelogram_jacobian(h, gamma, alpha_star, P_star, period)
% Given a parallelogram defined by 'h' and 'gamma', the function computes the Jacobian 
% at a certain point using central differencing for the derivatives.
% The point of interest (P_star, alpha_star) is part of a periodic orbit with period 'period'.

% Jacobian
% [df/d(alpha)        df/dP
%  dg/d(alpha)        dg/d

% where alpha_{n+1} = f(alpha_n, P_{n}) and P_{n+1} = g(alpha_n, P_{n})

eps = 1e-6; % For the finite difference. We use the central finite difference method

% How many iteations should we go ahead? (How many compositions of the
% maps?)
N = period;


% df/d(alpha) = (f(alpha + eps, P) - f(alpha - eps, P)) / (2*eps)
[s df_da_big_alpha df_da_big_P] = parallelogram_map(h, gamma, alpha_star+eps, P_star, N);       % f(alpha+eps, P)
[s df_da_small_alpha df_da_small_P] = parallelogram_map(h, gamma, alpha_star-eps, P_star, N);   % f(alpha-eps, P)
df_da = (df_da_big_alpha(N+1)-df_da_small_alpha(N+1))/(2*eps);               % approximately df/d(alpha)

% df/dP = (f(alpha, P + eps) - f(alpha, P - eps)) / (2*eps)
[s df_dP_big_alpha df_dP_big_P] = parallelogram_map(h, gamma, alpha_star, P_star+eps, N);       % f(alpha, P+eps)
[s df_dP_small_alpha df_dP_small_P] = parallelogram_map(h, gamma, alpha_star, P_star-eps, N);   % f(alpha, P-eps)
df_dP = (df_dP_big_alpha(N+1)-df_dP_small_alpha(N+1))/(2*eps);               % approximately df/dP

%  dg/d(alpha) = (g(alpha + eps, P) - g(alpha - eps, P)) / (2*eps)
[s dg_da_big_alpha dg_da_big_P] = parallelogram_map(h, gamma, alpha_star+eps, P_star, N);       % g(alpha+eps, P)
[s dg_da_small_alpha dg_da_small_P] = parallelogram_map(h, gamma, alpha_star-eps, P_star, N);   % g(alpha-eps, P)
dg_da = (dg_da_big_P(N+1)-dg_da_small_P(N+1))/(2*eps);                       % approximately dg/d(alpha)

%  dg/dP = (g(alpha, P + eps) - g(alpha, P - eps)) / (2*eps)
[s dg_dP_big_alpha dg_dP_big_P] = parallelogram_map(h, gamma, alpha_star, P_star+eps, N);       % g(alpha, P+eps)
[s dg_dP_small_alpha dg_dP_small_P] = parallelogram_map(h, gamma, alpha_star, P_star-eps, N);   % g(alpha, P-eps)
dg_dP = (dg_dP_big_P(N+1)-dg_dP_small_P(N+1))/(2*eps);                       % approximately dg/dP


J = [df_da df_dP; dg_da dg_dP];

% eig(J)

end

