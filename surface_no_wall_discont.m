function [] = surface_no_wall_discont(X, Y, Z, C)
%UNTITLED3 Summary of this function goes here
    h = surf(X,Y,Z, C);

    hp = patch(surf2patch(h), 'LineStyle', 'None');
    delete(h);
    V = get(hp,'Vertices');
    F = get(hp,'Faces');
    % Set the Alpha to be zero when the "slope" of a face is beyond a threshold
    A = ones(prod(size(Z)-1),1);
    thresh = 1.0; % increase to make more tolerant
    for n = 1:size(F,1)
        z = V(F(n,:),3);
        dz = max(max(abs(bsxfun(@minus,z,z'))));
        if dz > thresh;
            A(n) = 0;
        end
    end
    set(hp,'FaceVertexAlphaData',A,'edgealpha','flat');
    shading faceted;
    alpha flat


end

