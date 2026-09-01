function [vel_xyz] = convert_beam_to_xyz_velocity(beam1,beam2,beam3,beam4,beamAngle,excludeBeam)
    
    % matrix to convert from beam velocities to xyz velocity in instrument
    % coordinate frame
    if excludeBeam == 0
        conversion_matrix = [1/(2*sind(beamAngle)) 0 -1/(2*sind(beamAngle)) 0;
            0 -1/(2*sind(beamAngle)) 0 1/(2*sind(beamAngle));
            1/(4*cosd(beamAngle)) 1/(4*cosd(beamAngle)) 1/(4*cosd(beamAngle)) 1/(4*cosd(beamAngle))];

        vel_xyz = conversion_matrix*[beam1; beam2; beam3; beam4];

    elseif excludeBeam == 1
        conversion_matrix = [1/(2*sind(beamAngle)) -1/sind(beamAngle) 1/(2*sind(beamAngle));
            -1/(2*sind(beamAngle)) 0 1/(2*sind(beamAngle));
            1/(2*cosd(beamAngle)) 0 1/(2*cosd(beamAngle))];

        vel_xyz = conversion_matrix*[beam2; beam3; beam4];

    elseif excludeBeam == 2
        conversion_matrix = [1/(2*sind(beamAngle)) -1/(2*sind(beamAngle)) 0;
            -1/(2*sind(beamAngle)) -1/(2*sind(beamAngle)) 1/sind(beamAngle);
            1/(2*cosd(beamAngle)) 1/(2*cosd(beamAngle)) 0];

        vel_xyz = conversion_matrix*[beam1; beam3; beam4];

    elseif excludeBeam == 3
        conversion_matrix = [1/(sind(beamAngle)) -1/(2*sind(beamAngle)) -1/(2*sind(beamAngle));
            0 -1/(2*sind(beamAngle)) 1/(2*sind(beamAngle));
            0 1/(2*cosd(beamAngle)) 1/(2*cosd(beamAngle))];

        vel_xyz = conversion_matrix*[beam1; beam2; beam4];

    elseif excludeBeam == 4
        conversion_matrix = [1/(2*sind(beamAngle)) 0 -1/(2*sind(beamAngle));
            1/(2*sind(beamAngle)) -1/sind(beamAngle) 1/(2*sind(beamAngle));
            1/(2*cosd(beamAngle)) 0 1/(2*cosd(beamAngle))];

        vel_xyz = conversion_matrix*[beam1; beam2; beam3];
    end
end