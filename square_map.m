function [alpha, position] = square_map(alpha0, position0, N)


alpha = zeros(1, N+1); side = zeros(1, N+1); position = zeros(1, N+1);
% alpha (0, pi)        floor(P) {0,1,2,3}    P  [0,4)

%%% Initial Conditions
side(1) = floor(position0); position(1) = mod(position0,4); alpha(1) = alpha0;


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

end

