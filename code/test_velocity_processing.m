% specify location
latitude = 71.3274;
longitude = -156.8791;

% specify magnetic north declination
declination = 10.5; % degrees, positive if magnetic north is east of true north
% calculated for 2025-02-01 using https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml

% specify beam angle (degrees)
beam_angle = 25;

% processing parameters
mincorr = 40; % units are %, this is correlation of 0.4. Value taken from Jim's processing code

tic
%% load data
file_number = 41;
load(['C:/Users/hjohn/Documents/work/utqiagvik_mooring/Matlab_Format_with_coord_transforms/S106174A002_let_s_go_' num2str(file_number) '.mat']);

% set plotting limits
datenum_lims = [min(Data.Average_Time) max(Data.Average_Time)];

% initialize average counter (across multiple files)
acounter = 1;

% initialize structure
sigAverage = struct;

% specify depth grid (but this is actually height above the bottom!
% Confusing!
doff = 0.25; % Steve's estimate
depth_cell_size = 0.5;
depth_cells = (1:depth_cell_size:40);

% Find indices for each ensemble
firstindex = find(Data.Average_EnsembleCount==1);
spacing = mean(diff(firstindex));
na = length(firstindex);

% set velocities and backscatters to NaN where correlation is low
trim_ind = Data.Average_CorBeam1 < mincorr;
Data.Average_VelBeam1(trim_ind) = NaN;
Data.Average_AmpBeam1(trim_ind) = NaN;

trim_ind = Data.Average_CorBeam2 < mincorr;
Data.Average_VelBeam2(trim_ind) = NaN;
Data.Average_AmpBeam2(trim_ind) = NaN;

trim_ind = Data.Average_CorBeam3 < mincorr;
Data.Average_VelBeam3(trim_ind) = NaN;
Data.Average_AmpBeam3(trim_ind) = NaN;

trim_ind = Data.Average_CorBeam4 < mincorr;
Data.Average_VelBeam4(trim_ind) = NaN;
Data.Average_AmpBeam4(trim_ind) = NaN;

% Loop through averages and process
for a = 1:na
    AvgInd = firstindex(a) + (1:spacing) - 1;
    AvgInd( AvgInd > length(Data.Average_Time) ) = [];

    % meta data
    sigAverage(acounter).time = Data.Average_Time(AvgInd(1));
    sigAverage(acounter).lat = latitude;
    sigAverage(acounter).lon = longitude;
    sigAverage(acounter).watertemp = median(Data.Average_Temperature(AvgInd));
    sigAverage(acounter).depth = median(Data.Average_Pressure(AvgInd));
    sigAverage(acounter).z = depth_cells;
    
    % initialize arrays to store ping velocities
    [east,north,up,vel_x,vel_y,vel_z,backscatter1,backscatter2,backscatter3,backscatter4] = deal(nan(length(AvgInd),length(depth_cells)));

    % iterate over each individual ping
    for i = 1:length(AvgInd)

        % create rotation matrices
        R_heading = create_rotation_matrix('z',-Data.Average_Heading(AvgInd(i)));
        R_pitch = create_rotation_matrix('y',-Data.Average_Pitch(AvgInd(i)));
        R_roll = create_rotation_matrix('x',Data.Average_Roll(AvgInd(i)));
        R_enu = create_rotation_matrix('z',90 - declination);

        % rotation matrix to convert xyz in instrument coordinates to ENU
        R_xyz_to_enu = R_enu*R_heading*R_pitch*R_roll;

        % calculate vectors specifying beam positions in xyz coordinates
        % (instrument coordinate frame)
        beam1_xyz = [sind(beam_angle); 0; cosd(beam_angle)];
        beam2_xyz = [0; -sind(beam_angle); cosd(beam_angle)];
        beam3_xyz = [-sind(beam_angle); 0; cosd(beam_angle)];
        beam4_xyz = [0; sind(beam_angle); cosd(beam_angle)];

        % transform beam position vectors to ENU coordinates
        beam1_enu = R_xyz_to_enu*beam1_xyz;
        beam2_enu = R_xyz_to_enu*beam2_xyz;
        beam3_enu = R_xyz_to_enu*beam3_xyz;
        beam4_enu = R_xyz_to_enu*beam4_xyz;

        % create array of cell position along beam
        cell_distance_along_beam = (Config.Average_BlankingDistance + Config.Average_CellSize.*(1:double(Data.Average_NCells(AvgInd(i)))))./cosd(beam_angle);
        
        % calculate depth cell position along each individual beam
        beam1_depth_cell_positions_along_beam = (depth_cells-doff)./beam1_enu(3);
        beam2_depth_cell_positions_along_beam = (depth_cells-doff)./beam2_enu(3);
        beam3_depth_cell_positions_along_beam = (depth_cells-doff)./beam3_enu(3);
        beam4_depth_cell_positions_along_beam = (depth_cells-doff)./beam4_enu(3);
        
        % create beam velocity interpolants
        beam1_vel_interpolant = griddedInterpolant(cell_distance_along_beam,Data.Average_VelBeam1(AvgInd(i),:),'linear');
        beam2_vel_interpolant = griddedInterpolant(cell_distance_along_beam,Data.Average_VelBeam2(AvgInd(i),:),'linear');
        beam3_vel_interpolant = griddedInterpolant(cell_distance_along_beam,Data.Average_VelBeam3(AvgInd(i),:),'linear');
        beam4_vel_interpolant = griddedInterpolant(cell_distance_along_beam,Data.Average_VelBeam4(AvgInd(i),:),'linear');

        % calculate beam velocity at each depth
        beam1_vel_at_depth_cells = beam1_vel_interpolant(beam1_depth_cell_positions_along_beam);
        beam2_vel_at_depth_cells = beam2_vel_interpolant(beam2_depth_cell_positions_along_beam);
        beam3_vel_at_depth_cells = beam3_vel_interpolant(beam3_depth_cell_positions_along_beam);
        beam4_vel_at_depth_cells = beam4_vel_interpolant(beam4_depth_cell_positions_along_beam);

        % create beam amplitude interpolants
        beam1_amp_interpolant = griddedInterpolant(cell_distance_along_beam,Data.Average_AmpBeam1(AvgInd(i),:),'linear');
        beam2_amp_interpolant = griddedInterpolant(cell_distance_along_beam,Data.Average_AmpBeam2(AvgInd(i),:),'linear');
        beam3_amp_interpolant = griddedInterpolant(cell_distance_along_beam,Data.Average_AmpBeam3(AvgInd(i),:),'linear');
        beam4_amp_interpolant = griddedInterpolant(cell_distance_along_beam,Data.Average_AmpBeam4(AvgInd(i),:),'linear');

        % calculate beam amplitude at each depth and populate arrays
        backscatter1(i,:) = beam1_amp_interpolant(beam1_depth_cell_positions_along_beam);
        backscatter2(i,:) = beam2_amp_interpolant(beam2_depth_cell_positions_along_beam);
        backscatter3(i,:) = beam3_amp_interpolant(beam3_depth_cell_positions_along_beam);
        backscatter4(i,:) = beam4_amp_interpolant(beam4_depth_cell_positions_along_beam);

        % trim each beam to exclude sidelobe effects
        max_depth = Data.Average_Pressure(AvgInd(i)).*min([beam1_enu(3),beam3_enu(3),beam4_enu(3)]);
        trim_ind = depth_cells > (max_depth - depth_cell_size);
        beam1_vel_at_depth_cells(trim_ind) = NaN;
        beam3_vel_at_depth_cells(trim_ind) = NaN;
        beam4_vel_at_depth_cells(trim_ind) = NaN;

        max_depth = Data.Average_Pressure(AvgInd(i)).*beam1_enu(3) - depth_cell_size;
        trim_ind = depth_cells > max_depth;
        beam1_vel_at_depth_cells(trim_ind) = NaN;
        backscatter1(i,trim_ind) = NaN;

        max_depth = Data.Average_Pressure(AvgInd(i)).*beam2_enu(3) - depth_cell_size;
        trim_ind = depth_cells > max_depth;
        beam2_vel_at_depth_cells(trim_ind) = NaN;
        backscatter2(i,trim_ind) = NaN;

        max_depth = Data.Average_Pressure(AvgInd(i)).*beam3_enu(3) - depth_cell_size;
        trim_ind = depth_cells > max_depth;
        beam3_vel_at_depth_cells(trim_ind) = NaN;
        backscatter3(i,trim_ind) = NaN;

        max_depth = Data.Average_Pressure(AvgInd(i)).*beam4_enu(3) - depth_cell_size;
        trim_ind = depth_cells > max_depth;
        beam4_vel_at_depth_cells(trim_ind) = NaN;
        backscatter4(i,trim_ind) = NaN;

        % convert beam velocities to xyz velocity
        vector_velocity_xyz = convert_beam_to_xyz_velocity(beam1_vel_at_depth_cells,beam3_vel_at_depth_cells,beam4_vel_at_depth_cells,beam_angle);

        % convert velocity to ENU
        vector_velocity_enu = R_xyz_to_enu*vector_velocity_xyz;

        % extract velocity components and populate arrays
        east(i,:) = vector_velocity_enu(1,:);
        north(i,:) = vector_velocity_enu(2,:);
        up(i,:) = vector_velocity_enu(3,:);

    end

    % calculate average velocities over averaging interval
    sigAverage(acounter).east = mean(east,1,'omitnan');
    sigAverage(acounter).north = mean(north,1,'omitnan');
    sigAverage(acounter).up = mean(up,1,'omitnan');

    % calculate average backscatter
    sigAverage(acounter).backscatter1 = mean(backscatter1,1,'omitnan');
    sigAverage(acounter).backscatter2 = mean(backscatter2,1,'omitnan');
    sigAverage(acounter).backscatter3 = mean(backscatter3,1,'omitnan');
    sigAverage(acounter).backscatter4 = mean(backscatter4,1,'omitnan');

    % increment average counter
    acounter = acounter + 1;
end

toc

% extract data into useable matrices and vectors
[u,v,w,backscatter1,backscatter2,backscatter3,backscatter4] = deal(zeros(length(sigAverage),length(depth_cells)));
[time,watertemp] = deal(zeros(length(sigAverage),1));
for i = 1:length(sigAverage)
    time(i) = sigAverage(i).time;
    u(i,:) = sigAverage(i).east;
    v(i,:) = sigAverage(i).north;
    w(i,:) = sigAverage(i).up;
    watertemp(i) = sigAverage(i).watertemp;
    backscatter1(i,:) = sigAverage(i).backscatter1;
    backscatter2(i,:) = sigAverage(i).backscatter2;
    backscatter3(i,:) = sigAverage(i).backscatter3;
    backscatter4(i,:) = sigAverage(i).backscatter4;
end
%%
figure(1);
colormap(cmocean('balance'));

tiledlayout(3,1);
nexttile;
pcolor(time,depth_cells,u','EdgeColor','none');
cb = colorbar;
cb.Label.String = 'East (m/s)';
clim(max(abs(u),[],'all')*[-1 1]);
datetick('x');

nexttile;
pcolor(time,depth_cells,v','EdgeColor','none');
cb = colorbar;
cb.Label.String = 'North (m/s)';
clim(max(abs(v),[],'all')*[-1 1]);
datetick('x');

nexttile;
pcolor(time,depth_cells,w','EdgeColor','none');
cb = colorbar;
cb.Label.String = 'Up (m/s)';
clim(max(abs(w),[],'all')*[-1 1]);
datetick('x');

figure(2);
colormap(cmocean('amp'));
tiledlayout(4,1);

nexttile;
pcolor(time,depth_cells,backscatter1','edgecolor','none');
cb = colorbar;
cb.Label.String = 'Beam1 Amp';
datetick('x');

nexttile;
pcolor(time,depth_cells,backscatter2','edgecolor','none');
cb = colorbar;
cb.Label.String = 'Beam2 Amp';
datetick('x');

nexttile;
pcolor(time,depth_cells,backscatter3','edgecolor','none');
cb = colorbar;
cb.Label.String = 'Beam3 Amp';
datetick('x');

nexttile;
pcolor(time,depth_cells,backscatter4','edgecolor','none');
cb = colorbar;
cb.Label.String = 'Beam4 Amp';
datetick('x');