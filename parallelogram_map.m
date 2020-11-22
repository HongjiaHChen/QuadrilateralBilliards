% Map for parallelogram billiard

% The base length of the parallelogram is fixed at 1 and height is
% determined by the user input 'h'. The bottom left angle of the
% parallelogram is determined by 'gamma'. These two parameters uniquely
% define the parallelogram. User then inputs where the billiard ball
% initially starts at with 'alpha0' and 'position0' and how many iterations
% to compute with 'N'.

function [side alpha position] = parallelogram_map(h, gamma, alpha0, position0, N)


alpha = zeros(1, 1); side = zeros(1, 1); position = zeros(1, 1);

%%% Initial Conditions
side(1) = psuedo_floor(h, gamma, position0); position(1) = position0; alpha(1) = alpha0;

% It is not good to test for equality using == with floating point values,
% we need to introduce a thershold value, eps.
eps = 1e-13;


for i=1:N
    %disp(i)
    
    %position
    %side
    
    x_i = position(i)-side(i); % makes the code a bit more tidy

    % What side do we end up on? Then what angle and position?
    
    if (0 <= position(i) && position(i) <= 1) || (1+h/sin(gamma) <= position(i) && position(i) <= 2+h/sin(gamma))% Top or bottom side
        
        if 0 < alpha(i) && alpha(i) < atan(h/(1-x_i+h/tan(gamma)))  % right adjacent side
            %disp('Bottom/Top, fired right')
            
            side(i+1) = mod(side(i) + 1, 2*(1+h/sin(gamma)));
            
            if abs(side(i+1) - 2*(1+h/sin(gamma))) < eps
                side(i+1) = 0;                              % TRYING TO FIX OUR ERROR
            end
            
            position(i+1) = side(i+1) + sin(alpha(i))/sin(gamma-alpha(i)) * (1-x_i);
            alpha(i+1) = gamma - alpha(i);
            
        %elseif atan(h/(1-x_i+h/tan(gamma))) < alpha(i) && alpha(i) < (pi - atan(h/(h/tan(gamma) - x_i))) % opposite side
        elseif atan(h/(1-x_i+h/tan(gamma))) < alpha(i) && alpha(i) < (mod(atan(h/(h/tan(gamma) - x_i)),pi)) % opposite side    
            
            %disp('Bottom/Top, fired opp')
            
            side(i+1) = mod(side(i) + 1 + h/sin(gamma), 2*(1+h/sin(gamma)));
            
            if abs(side(i+1) - 2*(1+h/sin(gamma))) < eps
                side(i+1) = 0;                              % TRYING TO FIX OUR ERROR
            end
            
            position(i+1) = side(i+1) + 1 - x_i + h/tan(gamma) -h/tan(alpha(i));
            alpha(i+1) = pi - alpha(i);
            
            
        %elseif (pi - atan(h/(h/tan(gamma) - x_i))) < alpha(i) && alpha(i) < pi        % left adjacent side
        elseif (mod(atan(h/(h/tan(gamma) - x_i)),pi)) < alpha(i) && alpha(i) < pi        % left adjacent side
            %disp('Bottom/Top, fired left')
            
            side(i+1) = mod(side(i) + 2 + h/sin(gamma), 2*(1+h/sin(gamma)));
            
            if abs(side(i+1) - 2*(1+h/sin(gamma))) < eps
                side(i+1) = 0;                              % TRYING TO FIX OUR ERROR
            end
            
            position(i+1) = side(i+1) + h/sin(gamma) - sin(alpha(i))/sin(alpha(i)-gamma) * x_i;
            alpha(i+1) = pi - alpha(i) + gamma;
            
        else                                                    % We hit a corner
            warning('How did we get here?')
            break
            
        end
        
        
    elseif (1 <= position(i) && position(i) <= 1+h/sin(gamma)) || (2+h/sin(gamma) <= position(i) && position(i) <= 2+2*h/sin(gamma))% Angled left/right sides
        %disp('Angled')
        
        tan2_g = tan(gamma/2)^2;
        % vertex_exp = mod(2*atan(0.5*cot(gamma/2)*(-y*tan2_g+sqrt((tan2_g*(y+1)+y-1)^2+4*tan2_g)-tan2_g-y+1)),pi)
        vertex_exp1 = mod(2*atan(0.5*cot(gamma/2)*(-(h/sin(gamma)-x_i)*tan2_g+sqrt((tan2_g*((h/sin(gamma)-x_i)+1)+(h/sin(gamma)-x_i)-1)^2+4*tan2_g)-tan2_g-(h/sin(gamma)-x_i)+1)),pi);
        vertex_exp2 = mod(2*atan(0.5*cot(gamma/2)*(x_i*tan2_g+sqrt((tan2_g*(-x_i+1)-x_i-1)^2+4*tan2_g)-tan2_g+x_i+1)),pi);
        
        %if 0 < alpha(i) && alpha(i) < acot((cos(gamma)+h/sin(gamma)-x_i)/sin(gamma))  % right adjacent side
        if 0 < alpha(i) && alpha(i) < vertex_exp1  % right adjacent side
            %disp('Angled, fired right')
            
            side(i+1) = mod(side(i) + h/sin(gamma), 2*(1+h/sin(gamma)));
            
            if abs(side(i+1) - 2*(1+h/sin(gamma))) < eps
                side(i+1) = 0;                              % TRYING TO FIX OUR ERROR
            end
            
            position(i+1) = side(i+1) + sin(alpha(i))/sin(alpha(i)+gamma) * (h/sin(gamma)-x_i);
            alpha(i+1) = pi - alpha(i) - gamma;
            
        %elseif acot((cos(gamma)+h/sin(gamma)-x_i)/sin(gamma)) < alpha(i) && alpha(i) < (pi + acot((cos(gamma)-x_i)/sin(gamma))) % opposite side
        %elseif acot((cos(gamma)+h/sin(gamma)-x_i)/sin(gamma)) < alpha(i) && alpha(i) < mod(acot((cos(gamma)-x_i)/sin(gamma)),pi) % opposite side
        elseif vertex_exp1 < alpha(i) && alpha(i) < vertex_exp2 % opposite side
            %disp('Angled, fired opp')
            
            side(i+1) = mod(side(i) + 1 + h/sin(gamma), 2*(1+h/sin(gamma)));
            
            if abs(side(i+1) - 2*(1+h/sin(gamma))) < eps
                side(i+1) = 0;                              % TRYING TO FIX OUR ERROR
            end
            
            %position(i+1) = side(i+1) + h/sin(gamma) - x_i - sin(gamma+alpha(i))/sin(alpha(i));
            
            position(i+1) = side(i+1) + abs(h/sin(gamma) - x_i - sin(gamma+alpha(i))/sin(alpha(i)));
            alpha(i+1) = pi - alpha(i);
            
        %elseif (pi + acot((cos(gamma)-x_i)/sin(gamma))) < alpha(i) && alpha(i) < pi        % left adjacent side
        %elseif mod(acot((cos(gamma)-x_i)/sin(gamma)),pi) < alpha(i) && alpha(i) < pi        % left adjacent side
        elseif vertex_exp2 < alpha(i) && alpha(i) < pi        % left adjacent side
            
            %disp('Angled, fired left')
            
            side(i+1) = mod(side(i) + 1 + 2*h/sin(gamma), 2*(1+h/sin(gamma)));
                        
            if abs(side(i+1) - 2*(1+h/sin(gamma))) < eps
                side(i+1) = 0;                              % TRYING TO FIX OUR ERROR
            end
            
            position(i+1) = side(i+1) + 1 + x_i * sin(alpha(i))/sin(gamma+alpha(i));
            alpha(i+1) = 2*pi - gamma - alpha(i);
        
        else                                                    % We hit a corner
            warning('How did we get here?')
            break
        end
    else
        warning('Not in range')
            break
        %N = i-1; % for plotting purposes
        %side = side(1:N); position = position(1:N); alpha = alpha(1:N);
        break
    end
end
end


function side = psuedo_floor(h, gamma, raw_pos)
    %The usual floor does not work as we are not working with integers. 
    %The possible side indices are in {0, 1, 1 + h/sin(gamma), 2 + h/sin(gamma)}.
    %The raw position raw_pos lies in the range (0, 2+2 h/sin(gamma)).
    
    if 0 <= raw_pos && raw_pos <= 1
        side = 0;
    elseif 1 < raw_pos && raw_pos <= 1 + h/sin(gamma)
        side = 1;
    elseif 1 + h/sin(gamma) < raw_pos && raw_pos <= 2 + h/sin(gamma)
        side = 1 + h/sin(gamma);
    elseif 2 + h/sin(gamma) < raw_pos && raw_pos <= 2 + 2 * h/sin(gamma)
        side = 2 + h/sin(gamma);
    else
        disp(sprintf('raw_pos not in range (0, %0.2f)', 2 + 2 * h/sin(gamma)))
    end
end