function R = create_rotation_matrix(axis,angle)

    % Creates matrix describing rotation in cartesian coordinates about the
    % specified axis (x, y, or z) and by the specified angle (in degrees).
    % A positive angle corresponds to a rotation of a vector about the axis 
    % as indicated by the right-hand rule - for instance, R_z(angle) is
    % counter-clockwise for a positive angle. Alternately, this describes a
    % rotation of the coordinate system in the opposite direction.

    % Hayden Johnson, 2025-12-15
    
    if strcmp(axis,'x') | strcmp(axis,'X')
        R = [1 0 0;
            0 cosd(angle) -sind(angle);
            0 sind(angle) cosd(angle)];
    elseif strcmp(axis,'y') | strcmp(axis,'Y')
        R = [cosd(angle) 0 sind(angle);
            0 1 0;
            -sind(angle) 0 cosd(angle)];
    elseif strcmp(axis,'z') | strcmp(axis,'Z')
        R = [cosd(angle) -sind(angle) 0;
            sind(angle) cosd(angle) 0;
            0 0 1];
    end
end