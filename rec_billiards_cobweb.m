% MATLAB needs the first function to be the same name as the filename

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We draw a plot on the x_{n} - x_{n+1} plane.

% Our Difference Equation for Rectangular Billiards is on paper, it is too
% complicated to note here. Our rectangle has horizontal length L_H and
% vertical length L_V.

% alpha0 = inital angle ((0, pi) - angle at which the billiard leaves the 
                           % table with respect to right side of the table)

                
% position0 = intial position ((0, 2(L_H + L_V)) - the position at which the billiard 
                                % first touches the table)
% N = number of iterations

% All the above should be in vpa() form, i.e: symbolic form for extra
% accuracy. Float numbers are not enough.

% The types (symbols) of plots we can make have the iterates being:
% dots: 'o'
% connected lines: '-'
% dotted lines: 'o-'

% Use "o" to see the discontinuity and manual cobwebbing.
% Use "-" for easier visuals of periodicity.


function [side alpha position] = rec_billiards_cobweb(L_H, L_V, alpha0, position0, N, type, plot_flag)


alpha = zeros(1, N+1); side = zeros(1, N+1); position = zeros(1, N+1);
% alpha (0, pi)     floor(P) {0,L_H,L_H+L_V,2L_H+L_V}    P  [0,2(L_H + L_V))

%%% Initial Conditions
side(1) = psuedo_floor(position0, L_H, L_V);
position(1) = position0; alpha(1) = alpha0;

% It is not good to test for equality using == with floating point values,
% we need to introduce a thershold value, eps.
eps = 1e-13;

for i=1:N
    position(i) = vpa(position(i));
    side(i) = vpa(side(i));
    alpha(i) = vpa(alpha(i)); 
    
    x_i = vpa(position(i)-side(i)); % makes the code a bit more tidy
    
    if x_i < eps  % hit a vertex
        out_message = sprintf('Error: Last angle = %d radians. Last Position = %d.', alpha(i), position(i));
        disp(out_message)
        warning('We may have hit a corner and the trajectory has terminated')
        N = i; % for plotting purposes
        break
    end
    
    % What side do we end up on? Then what angle and position?
    
    % Bottom and top sides of length L_H
    if simplify(vpa(0) <= position(i) & position(i) < L_H) || simplify(L_H + L_V < position(i) & position(i) < 2 * L_H + L_V)
    
    %if abs(side(i)-0) < eps || abs(side(i) - (L_H + L_V)) < eps  % Bottom and top sides of length L_H

        if simplify(vpa(0) < alpha(i) & alpha(i) < atan(L_V/(L_H-x_i)))  % right adjacent side
            side(i+1) =  mod(vpa(side(i)) + L_H, L_H+L_H+L_V+L_V);
            if abs(side(i+1) - (L_H+L_H+L_V+L_V)) < eps
                side(i+1) = 0;                              % TRYING TO FIX OUR ERROR
            end
            alpha(i+1) = vpa(pi/2 - alpha(i));
            position(i+1) = side(i+1) + (L_H-x_i) * tan(vpa(alpha(i)));

        elseif simplify(atan(L_V/(L_H-x_i)) < alpha(i)) & simplify(alpha(i) < pi-atan(L_V/x_i)) % opposite side
            side(i+1) =  mod(vpa(side(i)) + L_H + L_V, L_H+L_H+L_V+L_V);
            if abs(side(i+1) - (L_H+L_H+L_V+L_V)) < eps
                side(i+1) = 0;                              % TRYING TO FIX OUR ERROR
            end
            alpha(i+1) = vpa(pi - alpha(i));
            position(i+1) = side(i+1) + L_H - x_i - L_V * cot(vpa(alpha(i)));

        elseif simplify(pi-atan(L_V/x_i) < alpha(i) & alpha(i) < vpa(pi))  % left adjacent side
            side(i+1) =  mod(vpa(side(i)) + L_H + L_H + L_V, L_H+L_H+L_V+L_V);
            if abs(side(i+1) - (L_H+L_H+L_V+L_V)) < eps
                side(i+1) = 0;                              % TRYING TO FIX OUR ERROR
            end
            alpha(i+1) = vpa(3 * pi/2 - alpha(i));
            position(i+1) = side(i+1) + L_V + x_i * tan(vpa(alpha(i)));
            
        else                                                    % We hit a corner
            warning('How did we get here?')
            break
        end
    % Left and right sides with length L_V
    elseif simplify(L_H < position(i) & position(i) < L_H + L_V) || simplify(2 * L_H + L_V < position(i) & position(i) < 2 * (L_H + L_V)) 
        
    %elseif abs(side(i)-L_H)<eps || abs(side(i)-(2*L_H+L_V))<eps  % Left and right sides with length L_V

        if simplify(vpa(0) < alpha(i) & alpha(i) < atan(L_H/(L_V-x_i)))  % right adjacent side
            side(i+1) =  mod(vpa(side(i)) + L_V, L_H+L_H+L_V+L_V);
            if abs(side(i+1) - (L_H+L_H+L_V+L_V)) < eps
                side(i+1) = 0;                              % TRYING TO FIX OUR ERROR
            end
            
            alpha(i+1) = vpa(pi/2 - alpha(i));
            position(i+1) = side(i+1) + (L_V-x_i) * tan(vpa(alpha(i)));

        elseif simplify(atan(L_H/(L_V-x_i)) < alpha(i)) && simplify(alpha(i) < pi-atan(L_H/x_i)) % opposite side
            side(i+1) =  mod(vpa(side(i)) + L_H + L_V, L_H+L_H+L_V+L_V);
            if abs(side(i+1) - (L_H+L_H+L_V+L_V)) < eps
                side(i+1) = 0;                              % TRYING TO FIX OUR ERROR
            end
            alpha(i+1) = vpa(pi - alpha(i));
            position(i+1) = side(i+1) + L_V - x_i - L_H * cot(vpa(alpha(i)));

        elseif simplify(pi-atan(L_H/x_i) < alpha(i) & alpha(i) < vpa(pi))  % left adjacent side
            side(i+1) =  mod(vpa(side(i)) + L_V + L_V + L_H, L_H+L_H+L_V+L_V);
            if abs(side(i+1) - (L_H+L_H+L_V+L_V)) < eps
                side(i+1) = 0;                              % TRYING TO FIX OUR ERROR
            end
            alpha(i+1) = vpa(3 * pi/2 - alpha(i));
            position(i+1) = side(i+1) + L_H + x_i * tan(vpa(alpha(i)));

        else                                                    % We hit a corner
            warning('How did we get here?')
            break
        end
    else
        N = i; % for plotting purposes
        out_message = sprintf('Error: Last angle = %d radians. Last Position = %d.', alpha(i), position(i));
        disp(out_message)
        disp('Error: We are not on a side!')
        break
    end
end

position = position(1:N); side = side(1:N); alpha = alpha(1:N); % delete the trailing zeros if necessary


if nargin <= 6 % user has not provided a plot_flag argument
    
    plot(position(1:(N-1)), position(2:(N)), type); hold on % depends on the symbol type the user wants
    plot([0:2*(L_H+L_V)], [0:2*(L_H+L_V)], 'r--')
    title(sprintf('Cobwebbing Rectangular Billiards with base length %0.2f and height %0.2f', L_H, L_V))
    xlabel('P_n')
    ylabel('P_{n+1}')
    
    xticks(double([0, L_H, L_H+L_V, 2*L_H + L_V, 2*(L_H + L_V)]))
    yticks(double([0, L_H, L_H+L_V, 2*L_H + L_V, 2*(L_H + L_V)]))
end

end


function [side] = psuedo_floor(raw_pos, L_H, L_V)
    % The usual floor does not work as L_H and L_V are not necessarily
    % integers and floor only gives integers. The possible side indices are
    % in {0, L_H, L_H + L_V, 2L_H + L_V}.
    % The raw position raw_pos lies in the range (0, 2L_H+2L_V).
    
    if 0 <= raw_pos && raw_pos <= L_H
        side = 0;
    elseif L_H < raw_pos && raw_pos <= L_H + L_V
        side = L_H;
    elseif L_H + L_V < raw_pos && raw_pos <= 2 * L_H + L_V
        side = L_H + L_V;
    elseif 2 * L_H + L_V < raw_pos && raw_pos <= 2 * (L_H + L_V)
        side = 2 * L_H + L_V;
    else
        disp(sprintf('raw_pos not in range (0, %0.2f)', 2*(L_H+L_V)))
    end
end
