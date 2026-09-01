function [vel_xyz] = convert_beam_to_xyz_velocity(beam1,beam2,beam3,beam4,beamAngle,excludeBeam)
    
    % Function to convert along-beam velocities to vector velocity in
    % instrument frame coordinates. Note that the beam velocities are
    % assumed to specified at consistent depths - i.e. instrument tilt
    % needs to be corrected for before calling this function.

    % Inputs:
    % beamN: (m/s) [1xn] Velocity measured by beam N. Positive is away from 
    % transducer. n is the number of depth bins.
    % beamAngle: (degrees) [scalar] Angle of the beams from the z axis in
    % the coordinate frame (25 degrees for Signature instruments).
    % excludeBeam: [int] Beam number to exclude for a 3-beam solution. If
    % 0, then a 4-beam solution is calculated.

    % Outputs:
    % vel_xyz: (m/s) [3xn] Vector velocity in the instrument coordinate
    % frame at each of the n depth bins of the input beam velocities.

    % Hayden Johnson, 2026-09-01

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