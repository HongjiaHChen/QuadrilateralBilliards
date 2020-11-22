% MATLAB needs the first function to be the same name as the filename

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We draw a plot on the x_{n} - x_{n+1} plane.

% Our Difference Equation for Square Billiards is on paper, it is too
% complicated to note here.

% alpha0 = inital angle ((0, pi) - angle at which the billiard leaves the 
                           % table with respect to right side of the table)

                
% position0 = intial position ((0, 4) - the position at which the billiard 
                                % first touches the table)
% N = number of iterations

% The types (symbols) of plots we can make have the iterates being:
% dots: 'o'
% connected lines: '-'
% dotted lines: 'o-'

% Use "o" to see the discontinuity and manual cobwebbing.
% Use "-" for easier visuals of periodicity.


function [side alpha position] = square_billiards_cobweb(alpha0, position0, N, type, plot_flag)


alpha = zeros(1, N+1); side = zeros(1, N+1); position = zeros(1, N+1);
% alpha (0, pi)        floor(P) {0,1,2,3}    P  [0,4)

%%% Initial Conditions
side(1) = floor(position0); position(1) = position0; alpha(1) = alpha0;


for i=1:N
    x_i = position(i)-side(i); % makes the code a bit more tidy
    
    % What side do we end up on? Then what angle and position?
    
    if 0 < alpha(i) && alpha(i) < acot(1-x_i)  % right adjacent side
        side(i+1) = mod(side(i) + 1, 4);
        alpha(i+1) = pi/2 - alpha(i);
        position(i+1) = side(i+1) + tan(alpha(i)) - x_i * tan(alpha(i));
        
    elseif acot(1-x_i) < alpha(i) && alpha(i) < pi-acot(x_i) % opposite side
        side(i+1) = mod(side(i) + 2, 4);
        alpha(i+1) = pi - alpha(i);
        position(i+1) = side(i+1) + 1 - x_i - cot(alpha(i));
        
    elseif pi-acot(x_i) < alpha(i) && alpha(i) < pi  % left adjacent side
        side(i+1) = mod(side(i) + 3, 4);
        alpha(i+1) = 3 * pi/2 - alpha(i);
        position(i+1) = side(i+1) + 1 + x_i * tan(alpha(i));
    else                                                    % We hit a corner
        % Which corner did we hit?
        if alpha(i) == 0
            position(i+1) = mod(side(i) + 1,4);
            side(i+1) = position(i+1);
            alpha(i+1) = 0;
        elseif alpha(i) == acot(1-x_i)
            position(i+1) = mod(side(i) + 2,4);
            side(i+1) = position(i+1);
            alpha(i+1) = 0;
        elseif alpha(i) == pi-acot(x_i)
            position(i+1) = mod(side(i) + 3,4);
            side(i+1) = position(i+1);
            alpha(i+1) = 0;
        elseif alpha(i) == pi
            position(i+1) = side(i);
            side(i+1) = position(i+1);
            alpha(i+1) = 0;
        end
        out_message = sprintf('Error: Last angle = %d radians. Last Position = %d.', alpha(i), position(i));
        disp(out_message)
        warning('We may have hit a corner and the trajectory has terminated')
        
        N = i+1; % for plotting purposes
        side = side(1:i+1); % delete the excess elements
        position = position(1:i+1); 
        alpha = alpha(1:i+1);
        break
    end

end

if nargin <= 4 % user has not provided a plot_flag argument
    plot(position(1:(N-1)), position(2:(N)), type), hold on % depends on the symbol type the user wants
    plot([0:4], [0:4], 'r--')
    title('Cobwebbing Square Billiards')
    xlabel('P_n')
    ylabel('P_{n+1}')

end

end

