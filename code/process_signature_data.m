% specify location
latitude = 71.3274;
longitude = -156.8791;

% specify magnetic north declination
declination = 10.5; % degrees, positive if magnetic north is east of true north
% value 10.5 calculated for 2025-02-01 using https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml

% specify beam angle (degrees) of transducer faces relative to instrument z axis
beam_angle = 25;

% processing parameters
mincorr = 40; % units are %, this is correlation of 0.4. Value taken from Jim's processing code
max_beam_tilt = 65; % maximum allowable angle of tilt from horizontal for an individual beam (degrees)
min_depth = 4; % meters, exclude data when instrument is shallower than this

% specify vertical bin grid
doff = 0.25; % Steve's estimate
vertical_bin_size = 0.5; % m
vertical_bins = (1:vertical_bin_size:40);

% specify destination directory for processed data
destination_dir = 'C:/Users/hjohn/Documents/work/utqiagvik_mooring/processed_data/';

% read file names
data_dir = 'C:/Users/hjohn/Documents/work/utqiagvik_mooring/Matlab_Format_with_coord_transforms/';
file_name_base = 'S106174A002_let_s_go_*.mat';

% END OF USER-SPECIFIED PARAMETERS
% -------------------------------------------------------------------------

% construct list of files to process
file_list = dir([data_dir file_name_base]);

% Caution: whether matlab understands "file9, file10" etc. correctly seems 
% to depend on the current lunar phase. The natsort package from the MATLAB
% file exchange fixes this problem.
file_list = natsortfiles(file_list);

% exclude average file (last in alphanumeric ordering)
file_list = file_list(1:end-1);

% initialize average counter (across multiple files)
acounter = 1;

% initialize sigAverage structure
sigAverage = struct;

for fi = 87:95 % length(file_list)

    % load file
    disp(['file ' num2str(fi) ' of ' num2str(length(file_list))])
    load([data_dir file_list(fi).name]);

    % create array of cell position along beam
    cell_distance_along_beam = (Config.Average_BlankingDistance + Config.Average_CellSize.*(1:double(Data.Average_NCells(1))))./cosd(beam_angle);
    
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
        sigAverage(acounter).z = vertical_bins;
        
        % initialize arrays to store ping velocities
        [east,north,up,backscatter1,backscatter2,backscatter3,backscatter4] = deal(nan(length(AvgInd),length(vertical_bins)));
    
        % iterate over each individual ping
        for i = 1:length(AvgInd)
    
            % create rotation matrices
            R_heading = create_rotation_matrix('z',-Data.Average_Heading(AvgInd(i)));
            R_pitch = create_rotation_matrix('y',-Data.Average_Pitch(AvgInd(i)));
            R_roll = create_rotation_matrix('x',Data.Average_Roll(AvgInd(i)));
            R_enu = create_rotation_matrix('z',90 - declination);
    
            % rotation matrix to convert xyz in instrument coordinates to ENU
            R_xyz_to_enu = R_enu*R_heading*R_pitch*R_roll;
    
            % initialize beam structure
            beam = struct;
    
            % create beam velocity interpolants
            beam(1).vel_interpolant = griddedInterpolant(cell_distance_along_beam,Data.Average_VelBeam1(AvgInd(i),:),'linear');
            beam(2).vel_interpolant = griddedInterpolant(cell_distance_along_beam,Data.Average_VelBeam2(AvgInd(i),:),'linear');
            beam(3).vel_interpolant = griddedInterpolant(cell_distance_along_beam,Data.Average_VelBeam3(AvgInd(i),:),'linear');
            beam(4).vel_interpolant = griddedInterpolant(cell_distance_along_beam,Data.Average_VelBeam4(AvgInd(i),:),'linear');
    
            % create beam amplitude interpolants
            beam(1).amp_interpolant = griddedInterpolant(cell_distance_along_beam,Data.Average_AmpBeam1(AvgInd(i),:),'linear');
            beam(2).amp_interpolant = griddedInterpolant(cell_distance_along_beam,Data.Average_AmpBeam2(AvgInd(i),:),'linear');
            beam(3).amp_interpolant = griddedInterpolant(cell_distance_along_beam,Data.Average_AmpBeam3(AvgInd(i),:),'linear');
            beam(4).amp_interpolant = griddedInterpolant(cell_distance_along_beam,Data.Average_AmpBeam4(AvgInd(i),:),'linear');
    
            % calculate vectors specifying beam positions in xyz coordinates
            % (instrument coordinate frame)
            beam(1).unit_vector_xyz = [sind(beam_angle); 0; cosd(beam_angle)];
            beam(2).unit_vector_xyz = [0; -sind(beam_angle); cosd(beam_angle)];
            beam(3).unit_vector_xyz = [-sind(beam_angle); 0; cosd(beam_angle)];
            beam(4).unit_vector_xyz = [0; sind(beam_angle); cosd(beam_angle)];
    
            for j = 1:length(beam)
                % transform beam position vectors to ENU coordinates
                beam(j).unit_vector_enu = R_xyz_to_enu*beam(j).unit_vector_xyz;
                
                % calculate vertical bin position along each individual beam
                beam(j).vertical_bin_positions_along_beam = (vertical_bins-doff)./beam(j).unit_vector_enu(3);
    
                % calculate beam velocity at each vertical bin
                beam(j).vel_at_vertical_bins = beam(j).vel_interpolant(beam(j).vertical_bin_positions_along_beam);
            
                % calculate beam amplitude (backscatter) at each vertical
                % bin
                beam(j).backscatter_at_vertical_bins = beam(j).amp_interpolant(beam(j).vertical_bin_positions_along_beam);
            
                % trim each beam to exclude sidelobe effects
                beam(j).max_height = Data.Average_Pressure(AvgInd(i)).*beam(j).unit_vector_enu(3) - vertical_bin_size;
                beam(j).min_height = doff + ((Config.Average_BlankingDistance + Config.Average_CellSize)./cosd(beam_angle)).*beam(j).unit_vector_enu(3) + vertical_bin_size;
                trim_ind = vertical_bins > beam(j).max_height | vertical_bins < beam(j).min_height;
                beam(j).vel_at_vertical_bins(trim_ind) = NaN;
                beam(j).backscatter_at_vertical_bins(trim_ind) = NaN;
            
                % check if any beams should be excluded based on tilt
                if beam(j).unit_vector_enu(3) < cosd(max_beam_tilt)
                    beam(j).exclude = true;
                else 
                    beam(j).exclude = false;
                end
            end
            
            exclude_array = [beam.exclude];
            if sum(exclude_array) > 1
                vector_velocity_xyz = NaN(3,length(vertical_bins));
            else
                if sum(exclude_array) == 0
                    exclude_beam = 0;
                elseif sum(exclude_array) == 1 
                    exclude_beam = find(exclude_array);
                end
                % convert beam velocities to instrument frame (xyz) velocity
                vector_velocity_xyz = convert_beam_to_xyz_velocity(beam(1).vel_at_vertical_bins,beam(2).vel_at_vertical_bins,beam(3).vel_at_vertical_bins,beam(4).vel_at_vertical_bins,beam_angle,exclude_beam);
            end
    
            % convert instrument frame (xyz) velocity to ENU
            vector_velocity_enu = R_xyz_to_enu*vector_velocity_xyz;
    
            % extract velocity components and populate arrays
            east(i,:) = vector_velocity_enu(1,:);
            north(i,:) = vector_velocity_enu(2,:);
            up(i,:) = vector_velocity_enu(3,:);
    
            % extract backscatter from each beam 
            backscatter1(i,:) = beam(1).backscatter_at_vertical_bins;
            backscatter2(i,:) = beam(2).backscatter_at_vertical_bins;
            backscatter3(i,:) = beam(3).backscatter_at_vertical_bins;
            backscatter4(i,:) = beam(4).backscatter_at_vertical_bins;
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
end

% remove averages where instrument was out of water
badavg = false(size(sigAverage));
for i = 1:length(sigAverage)
    if sigAverage(i).depth < min_depth
        badavg(i) = true;
    end
end
sigAverage(badavg) = [];

% extract data into useable matrices and vectors
[u,v,w,backscatter1,backscatter2,backscatter3,backscatter4] = deal(zeros(length(sigAverage),length(vertical_bins)));
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

%% create plot
figure(1);

% set velocity colour limits
combined_vel_array = [u v];
max_color = prctile(abs(combined_vel_array),99,'all');
vel_lims = max_color*[-1 1];

tiledlayout(3,1);
ax1 = nexttile;
pcolor(time,vertical_bins,u','EdgeColor','none');
colormap(ax1,cmocean('balance'));
cb = colorbar;
cb.Label.String = 'East (m/s)';
clim(vel_lims);
datetick('x');
ylabel ('z (m)');

ax2 = nexttile;
pcolor(time,vertical_bins,v','EdgeColor','none');
colormap(ax2,cmocean('balance'));
cb = colorbar;
cb.Label.String = 'North (m/s)';
clim(vel_lims);
datetick('x');
ylabel ('z (m)');

ax3 = nexttile;
pcolor(time,vertical_bins,backscatter4','edgecolor','none');
colormap(ax3,cmocean('deep'));
cb = colorbar;
cb.Label.String = 'Backscatter';
datetick('x');
ylabel ('z (m)');

%{
nexttile;
pcolor(time,vertical_bins,w','EdgeColor','none');
cb = colorbar;
cb.Label.String = 'Up (m/s)';
clim(max(abs(w),[],'all')*[-1 1]);
datetick('x');
ylabel ('z (m)');
%}