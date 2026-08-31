% specify ADCP geometry
beam_angle = 25;

% specify ADCP orientation
heading = 260;
pitch = 0.24;
roll = 52;
declination = 10.5;

water_depth = 42;

% create rotation matrices
R_heading = create_rotation_matrix('z',-heading);
R_pitch = create_rotation_matrix('y',-pitch);
R_roll = create_rotation_matrix('x',roll);

R_enu = create_rotation_matrix('z',90 - declination);

% This is the real transformation
R_combined = R_enu*R_heading*R_pitch*R_roll;

% wrong transformations for testing and verifying incorrectness
% R_combined = R_enu*R_roll*R_pitch*R_heading;
% R_combined = R_heading*R_pitch*R_roll*R_enu;

% calculate beam locations in ADCP coordinate system
z = (0:70);
x1 = z.*tand(beam_angle);
x3 = -z.*tand(beam_angle);
y2 = -z.*tand(beam_angle);
y4 = z.*tand(beam_angle);

[x2,x4,x5,y1,y3,y5] = deal(zeros(size(z)));

% create position vectors for beams
r1 = [x1; y1; z];
r2 = [x2; y2; z];
r3 = [x3; y3; z];
r4 = [x4; y4; z];
r5 = [x5; y5; z];

% calculate beam positions in rotated coordinate frame
r1_geo = R_combined*r1 - [0; 0; water_depth];
r2_geo = R_combined*r2 - [0; 0; water_depth];
r3_geo = R_combined*r3 - [0; 0; water_depth];
r4_geo = R_combined*r4 - [0; 0; water_depth];
r5_geo = R_combined*r5 - [0; 0; water_depth];

% plot beams in ADCP coordinates
figure(4);
plot3(x1,y1,z,'displayname','Beam 1');
hold on;
plot3(x2,y2,z,'displayname','Beam 2');
plot3(x3,y3,z,'displayname','Beam 3');
plot3(x4,y4,z,'displayname','Beam 4');
plot3(x5,y5,z,'color','black','displayname','Beam 5');
hold off;
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
legend();

% plot beams in ENU coordinates
figure(5);
plot3(r1_geo(1,:),r1_geo(2,:),r1_geo(3,:),'displayname','Beam 1');
hold on;
plot3(r2_geo(1,:),r2_geo(2,:),r2_geo(3,:),'displayname','Beam 2');
plot3(r3_geo(1,:),r3_geo(2,:),r3_geo(3,:),'displayname','Beam 3');
plot3(r4_geo(1,:),r4_geo(2,:),r4_geo(3,:),'displayname','Beam 4');
plot3(r5_geo(1,:),r5_geo(2,:),r5_geo(3,:),'color','black','displayname','Beam 5');
hold off;
grid on;
axis equal;
xlabel('east');
ylabel('north');
zlabel('depth');
legend();

