function [vel_xyz] = convert_beam_to_xyz_velocity(beam1,beam3,beam4,beamAngle)
    
    % matrix to convert from beam velocities to xyz velocity in instrument
    % coordinate frame
    conversion_matrix = [1/(2*sind(beamAngle)) -1/(2*sind(beamAngle)) 0;
        -1/(2*sind(beamAngle)) -1/(2*sind(beamAngle)) 1/sind(beamAngle);
        1/(2*cosd(beamAngle)) 1/(2*cosd(beamAngle)) 0];
    
    vel_xyz = conversion_matrix*[beam1; beam3; beam4];

end