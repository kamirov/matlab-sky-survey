function [t, ra, dec, fig] = sky_survey(el_min, t_test, fov, lat, lon, tuning, fig_idx, show_elevs)
% Generates a set of (RA, Dec) coordinates that encompass a region of the
% sky
% 
% Outputs:
%     - t               Trajectory datetimes (e.g. 01-Jan-2016 20:10:05)
%     - ra              Trajectory right ascensions [deg]
%     - dec             Trajectory declinations [deg]
%     - fig             Survey plot
%     
% Inputs:
%     - el_min          Minimum elevation [deg]
%     - t_test          Test time [min]
%     - fov             Field of view (half-cone) [deg]
%     - lat             Latitude [deg]
%     - lon             Longitude [deg]
%     - tuning          Tuning structure (optional)
%     - show_elevs      Boolean. If true, shows elevation profiles
%     - fig_idx         Figure index for main plot
% 
% Glitches:
%   - Doesn't work for lat = 0
%   - Test needs to be at least 1 timestep in length
% 
% Notes:
%   - Non-integer times are floored

tic


%% Quick Error Check

if lat == 0
    error('Enter a nonzero latitude.');
elseif t_test == 0
    error('Enter a positive nonzero test time.');
end


%% Initialize

fprintf('Surveying...');

n_sp = 100;             % # of sphere mesh points
t_test = floor(t_test);


% Radify
el_min = deg2rad(el_min);
lat = deg2rad(lat);
lon = deg2rad(lon);

% Setup plot
fig = figure(fig_idx);

clf
subplot(2,2,[1 3]);
axis equal;
alim = 1.1; 
axis([-alim alim -alim alim])
hold on;


%% Alt-az coordinates

% Celestial sphere
[A_sp, a_sp] = meshgrid(linspace(0,2*pi,n_sp), linspace(0,2*pi,n_sp));

% ROI separated into 2 semicircles (used later on)
[A, a] = meshgrid(linspace(pi, 2*pi, n_sp),linspace(el_min, pi/2, n_sp));
[A_e, a_e] = meshgrid(linspace(0, pi, n_sp),linspace(el_min, pi/2, n_sp));

% % Plot ROI semicircles
% figure(100);
% subplot(1,2,1);
% [x,y,z] = sph2cart(A, a, 1);
% surf(x,y,z);
% axis equal;
% axis([-alim alim -alim alim]);
% view([60 20]);
% subplot(1,2,2);
% [x,y,z] = sph2cart(A_e, a_e, 1);
% surf(x,y,z);
% axis equal;
% axis([-alim alim -alim alim]);
% view([60 20]);
% 
% return;

%% Convert alt-az to ra-dec

% ROI end-semicircle ra-dec is found at the test-end time
start_time = now;
utc = local_time_to_utc(start_time);
utc_e = local_time_to_utc(datenum(start_time + minutes(t_test)));
[r, d] = azel2radec(A, a, lat, lon, utc);
[r_e, d_e] = azel2radec(A_e, a_e, lat, lon, utc_e);

% Min/max are the outermost points
[r_min, ~] = azel2radec(3*pi/2, el_min, lat, lon, utc);
[r_max, ~] = azel2radec(pi/2, el_min, lat, lon, utc_e);
d_min = min(d(:));
d_max = max(d(:));

% Start-max / end-min are the innermost RAs of the start and end
% semicircles (we'll use this to get the region between the start and end
% semicircles)
if lat > 0
    ang = pi;
else
    ang = 0;
end
[r_start_max, ~] = azel2radec(ang, el_min, lat, lon, utc);
[r_end_min, ~] = azel2radec(ang, el_min, lat, lon, utc_e);

% Bounds fixes
if r_min < 0
    r_min = r_min + 2*pi;
    r_end_min = r_end_min + 2*pi;
    r_start_max = r_start_max + 2*pi;
end
if r_end_min < r_start_max
    r_end_min = r_end_min + 2*pi;
end
if r_max < r_min
    r_max = r_max + 2*pi;
end

% Intermediate region (between the start and end semicircles)
[r_int, d_int] = meshgrid(linspace(r_start_max, r_end_min), linspace(d_min, d_max, n_sp));


%% Plot ROI on celestial sphere

fig = figure(fig_idx);

% Celestial sphere
[x_sp, y_sp, z_sp] = sph2cart(A_sp, a_sp, 1);
surf(x_sp, y_sp, z_sp, 'FaceAlpha', 0.2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2);
% view(190, 10)
view([145 20])
% Plot start, intermediary, and end swath. Note, we can't merge them
% (the surface plot looks weird if we do; there's probably a fix)
[x_roi, y_roi, z_roi] = sph2cart(r, d, 1);          % Start
surf(x_roi, y_roi, z_roi, 'EdgeColor', 'none');
[x_roi, y_roi, z_roi] = sph2cart(r_int, d_int, 1);  % Intermediate
surf(x_roi, y_roi, z_roi, 'EdgeColor', 'none');
[x_roi, y_roi, z_roi] = sph2cart(r_e, d_e, 1);      % End
surf(x_roi, y_roi, z_roi, 'EdgeColor', 'none');

% Cartesian axes
line([0 1], [0 0], [0 0])
line([0 0], [0 1], [0 0])
line([0 0], [0 0], [0 1])


%% Add points

% # of points based on flattening the ROI, then seeing how many circles
% sized by the FOV will fit
num_pts = floor((180/(fov*2)) * 360/(fov*2));
pts = saff_sphere(num_pts);        % Uniformly distribute points [dec, ra]
pts(:,1) = pts(:,1) - pi/2;        % Declination correction

% Plot all points
% [x,y,z] = sph2cart(pts(:,2), pts(:,1), 1.05);
% scatter3(x,y,z, 'filled', 'k');
% return;


%% Remove points outside of ROI

pts_roi = [];   % [dec, ra]

% Increment our ROI in minutes over our test time to only select generated
% points that will be visible (minutes seems to be a good balance between
% performance and accuracy)
for i = 0:t_test
    utc_cur = datestr(datenum(utc, 'dd-mmm-yyyy HH:MM:SS') + minutes(i));
    [~, el_p] = radec2azel(pts(:,2), pts(:,1), lat, lon, utc_cur);

    for j = num_pts:-1:1
        num_pts_roi = size(pts_roi,1);
        if el_p(j) > el_min    
            % Add the point if either there are no other points or if this
            % point hasn't already been added (note, RA for any 2 points
            % will never be the same)
            if num_pts_roi <= 1 || ~any(pts_roi(:,2) == pts(j,2))
                pts_roi = [pts_roi ; pts(j,:)];
            end
        end
    end
end

% Bounds fix. There's a vectorized way to do this.
for i = 1:num_pts_roi
    if pts_roi(i,2) < r_min
        pts_roi(i,2) = pts_roi(i,2) + 2*pi;
    end
end

% Plot all ROI points
[x,y,z] = sph2cart(pts_roi(:,2), pts_roi(:,1), 1.05);
scatter3(x,y,z, 'filled', 'k');
% return

%% Get point elevation trajectories

elevs = [];
azs = [];
for i = 0:t_test
    utc_cur = datestr(datenum(utc, 'dd-mmm-yyyy HH:MM:SS') + minutes(i));
    [A_p, el_p] = radec2azel(pts_roi(:,2), pts_roi(:,1), lat, lon, utc_cur);

    for j = 1:num_pts_roi
        elevs(j,i+1) = el_p(j);
        azs(j,i+1) = A_p(j);
    end
end

% Max altitudes of each point and time of max altitudes
[max_alts, max_alt_ts] = max(elevs');
max_alt_ts = max_alt_ts - 1;  % MATLAB's 1-indexing offset

% 3rd point is the index (used for when we manipulate this array)
traj_max = [max_alt_ts', max_alts', (1:num_pts_roi)'];
traj_max = sortrows(traj_max, -2);  % Sort desc by elevation

% Plot elevation trajectories (only some points)
% figure(20)
% points = [1 10 20 30 40 25 27 35 36];
% scatter(max_alt_ts(points)', rad2deg(max_alts(points)),'k', 'filled');
% hold on;
% plot(0:t_test, rad2deg(elevs(points, :)))
% title('Point Elevation Trajectories')
% ylabel('Elevation [deg]')
% xlabel('Time [min]')
% legend('Peak elevations');
% return;

% figure(20)
% line([0, t_test], rad2deg([el_min, el_min]), 'LineStyle', ':', 'Color', 'k')
% hold on;
% scatter(max_alt_ts(:)', rad2deg(max_alts(:)),'k', 'filled');
% plot(0:t_test, rad2deg(elevs(:, :)))
% title('Point Elevation Trajectories')
% ylabel('Elevation [deg]')
% xlabel('Time [min]')
% legend('Minimum elevation', 'Peak elevations', 'Location', 'South');
% return;

%% Generate horizontal-coords trajectory

% clc

% Uniformly distribute time. Floor to allow us to get the corresponding
% elevation from the already found elevations matrix
t_uni = floor((0:num_pts_roi-1) * (t_test/(num_pts_roi-1)));
t_uni_tmp = t_uni;

% Pass 1 - Distribution
atraj_dist = [];     % [time, azimuth, elevation, pt_idx]
for i = 1:num_pts_roi
    % Find nearest uniform-time point
    [t_delta, t_idx] = min(abs(t_uni_tmp-traj_max(i,1)));
    az = azs(traj_max(i,3), t_uni_tmp(t_idx)+1);
    elev = elevs(traj_max(i,3), t_uni_tmp(t_idx)+1);
    
        atraj_dist = [atraj_dist
                      t_uni_tmp(t_idx), az, elev, traj_max(i,3)];

        t_uni_tmp(t_idx) = [];
end

% figure(20)
% line([0, t_test], rad2deg([el_min, el_min]), 'LineStyle', ':', 'Color', 'k')
% hold on;
% scatter(atraj_dist(:,1)', rad2deg(atraj_dist(:,3)),'r', 'filled');
% plot(0:t_test, rad2deg(elevs(:, :)))
% title('Point Elevation Trajectories')
% ylabel('Elevation [deg]')
% xlabel('Time [min]')
% legend('Minimum elevation', 'Elevation during imaging', 'Location', 'South');
% return;

if ~tuning.prioritize_high_elevations
    % Pass 2 - Fix low-elevations points (where possible)
    atraj_fix = atraj_dist;
    for i = 1:num_pts_roi
        if atraj_fix(i,3) < el_min
            pt1_old_t = atraj_fix(i,1);
            pt1_idx = atraj_fix(i,4);

            for j = 1:size(atraj_fix, 1)
                t = atraj_fix(j,1);
                pt1_new_el = elevs(pt1_idx, t+1);

                if pt1_new_el > el_min
                    pt1_new_t = t;
                    pt2_idx = atraj_fix(j,4);
                    pt2_new_el = elevs(pt2_idx, pt1_old_t+1);

                    if pt2_new_el > el_min
                        pt1_new = [pt1_new_t, azs(pt1_idx, pt1_new_t+1), ...
                                   elevs(pt1_idx, pt1_new_t+1), pt1_idx];
                        pt2_new = [pt1_old_t, azs(pt2_idx, pt1_old_t+1), ...
                                   elevs(pt2_idx, pt1_old_t+1), pt2_idx];

                        atraj_fix(i,:) = pt2_new;
                        atraj_fix(j,:) = pt1_new;

                        break;
                    end
                end
            end
        end
    end
    
    atraj = atraj_fix; % [time, az, elev, pt_idx]
else
    atraj = atraj_dist; % [time, az, elev, pt_idx]
end
   
% figure(20)
% line([0, t_test], rad2deg([el_min, el_min]), 'LineStyle', ':', 'Color', 'k')
% hold on;
% scatter(atraj(:,1)', rad2deg(atraj(:,3)),'r', 'filled');
% plot(0:t_test, rad2deg(elevs(:, :)))
% title('Point Elevation Trajectories')
% ylabel('Elevation [deg]')
% xlabel('Time [min]')
% legend('Minimum elevation', 'Elevation during imaging', 'Location', 'South');
% hold off
% return;

% Pass 3 - Remove low-elevation incurables
for i = num_pts_roi:-1:1
    if atraj(i, 3) < el_min
        atraj(i,:) = [];
    end
end
atraj = sortrows(atraj, 1);

% return

% figure(20)
% line([0, t_test], rad2deg([el_min, el_min]), 'LineStyle', ':', 'Color', 'k')
% hold on;
% scatter(atraj(:,1)', rad2deg(atraj(:,3)),'r', 'filled');
% plot(0:t_test, rad2deg(elevs(:, :)))
% title('Point Elevation Trajectories')
% ylabel('Elevation [deg]')
% xlabel('Time [min]')
% legend('Minimum elevation', 'Elevation during imaging', 'Location', 'South');
% return;


%% Adjust times to favor higher elevations

if tuning.timing.adjust
    atraj_adj = atraj;
    num_pts_atraj = size(atraj, 1);

    % Look at point to the left. If it's at a lower elevation, then shift
    % current point to the left to minimize lower elevation point's time
    for i = 2:num_pts_atraj
        cur_pt = atraj_adj(i,:);
        prev_pt = atraj_adj(i-1,:);

        if cur_pt(3) > prev_pt(3) && (cur_pt(1) - prev_pt(1) > tuning.timing.min_step)
            t_delta = max(tuning.timing.min_step, floor((cur_pt(1) - prev_pt(1)) / 2));
            atraj_adj(i,1) = cur_pt(1) - t_delta;
        end
    end

    % Look at point to the right. If it's at a higher elevation, then shift
    % current point to the right to minimize lower elevation point's time
    for i = 1:num_pts_atraj-1
        cur_pt = atraj_adj(i,:);
        next_pt = atraj_adj(i+1,:);

        if cur_pt(3) < next_pt(3) && (next_pt(1) - cur_pt(1) > tuning.timing.min_step)
            t_delta = max(tuning.timing.min_step, floor((next_pt(1) - cur_pt(1))/2));
            atraj_adj(i,1) = cur_pt(1) + t_delta;
        end
    end
    
    atraj = atraj_adj;
end


%% Get equatorial-coords trajectory

traj = [];  % [time, dec, ra]
num_pts_atraj = size(atraj, 1);
for i = num_pts_atraj : -1 : 1
    traj = [traj
            atraj(i,1), pts_roi(atraj(i,4), :)];
end

traj = sortrows(traj, 1);
num_pts_traj = num_pts_atraj;


if tuning.proximity.adjust
    %% Generate cost vector

    % Elevation matrix and updated azel trajectory
    els = [];
    atraj = [];     % [t, az, el]
    azs = [];
    k = 1;
    for t = traj(:,1)'
        utc_cur = datestr(datenum(utc, 'dd-mmm-yyyy HH:MM:SS') + minutes(t));
        [az_all_t, el_all_t] = radec2azel(traj(:,3), traj(:,2), lat, lon, utc_cur);
        
        els = [els
               el_all_t'];
           
        [az_t, el_t] = radec2azel(traj(k,3), traj(k,2), lat, lon, utc_cur);
        atraj = [atraj
                 t, az_t, el_t];
        k = k+1;
    end

    % Cost vector
    el_maxes = max(els);
    costs = [];
    for k = 1:num_pts_traj
        el_max = el_maxes(k);

        if atraj(k,3) < el_min
            c = 1;
        else
            c = 1 - (atraj(k,3)/el_max);
        end

        costs = [costs c];
    end


    %% Get path windows
    
    windows = cell(1,1);
    pt_idx = 1;
    window_idx = 1;
    while pt_idx < num_pts_traj + 1
        cur_sz = 1;
        cur_window = [];

        while cur_sz < tuning.proximity.window_max_size + 1
            if cur_sz == 1
                cur_pt = traj(pt_idx,:);
            else
                if costs(pt_idx) > tuning.proximity.window_cost_threshold ...
                   || cur_sz > tuning.proximity.window_max_size
                    break;
                else
                    cur_pt = traj(pt_idx,:);
                end
            end

            cur_window = [cur_window
                          cur_pt];
            cur_sz = cur_sz + 1;
            pt_idx = pt_idx + 1;
            
            if pt_idx > num_pts_traj
                break;
            end
        end

        windows{window_idx}= cur_window;
        window_idx = window_idx + 1;
    end

    % Plot windows
    for k = 1:numel(windows)
        window = windows{k};
        [x,y,z] = sph2cart(window(:,3), window(:,2), 1.05);
        scatter3(x,y,z, 'filled');  
    end
    
    %% Remake trajectory based on weighing ra-dec in each window

    traj_prox = [];
    num_windows = numel(windows);
    for k = 1:num_windows
        num_pts_window = size(windows{k}, 1);

        window_traj = [];
        for n = 1:num_pts_window
            if ~exist('cur_pt', 'var')
                window_traj(1,:) = windows{k}(1,:);
                cur_pt = window_traj(1,:);
                windows{k}(1,:) = [];
            else            
                % Get a proximity score for all points (weighed by RA
                % and dec factors)
                prox_costs = tuning.proximity.ra_weight*abs(cur_pt(3)-windows{k}(:,3)) ...
                             + tuning.proximity.dec_weight*abs(cur_pt(2)-windows{k}(:,2));

                nearest_idx = find(prox_costs == min(prox_costs), 1);
                cur_pt = windows{k}(nearest_idx,:);

                window_traj(n,:) = windows{k}(nearest_idx,:);
                windows{k}(nearest_idx,:) = [];
            end
        end

        traj_prox = [traj_prox
                      window_traj];
    end

    % At the end of this, the times are jumbled up, so sort them (without
    % sorting the whole trajectory)
    traj_radec_tmp = sortrows(traj_prox, 1);
    traj_prox(:,1) = traj_radec_tmp(:,1);

    traj = traj_prox;    
end


%% Plot point elevation trajectories

if show_elevs
    figure(fig_idx + 10000) % Arbitrary figure
    clf
    plot(0:t_test, rad2deg(elevs))
    title('Point Elevation Trajectories');
    xlabel('Time [min]');
    ylabel('Elevation [deg]');
    hold on
    line([0, t_test], rad2deg([el_min, el_min]), 'LineStyle', ':', 'Color', 'k')
    scatter(traj_max(:,1), rad2deg(traj_max(:,2)), 'filled', 'k')
    scatter(atraj_dist(:,1), rad2deg(atraj_dist(:,3)), 'filled', 'r')
    xlim([0 t_test])
    hold off
end


%% Plot path
fig = figure(fig_idx);

subplot(2,2, [1 3])


% [x,y,z] = sph2cart(traj(:,3), traj(:,2), 1.05);
% scatter3(x,y,z, 'filled', 'r');    

[x_path, y_path, z_path] = sph2cart(traj(:,3), traj(:,2), 1.05);
plot3(x_path, y_path, z_path, 'b')

plot_title = 'Sky Survey';

title(plot_title, 'FontSize', 12);

% Cosmetics
ax1 = gca;
yruler = ax1.YRuler;
yruler.Axle.Visible = 'off';
xruler = ax1.XRuler;
xruler.Axle.Visible = 'off';    
zruler = ax1.ZRuler;
zruler.Axle.Visible = 'off';    
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
set(gca,'ZTickLabel',[])
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'ZTick',[])
set(gca,'color','none')


%% Remove sub min-elevation points and plot

atraj_fine = [];    % [time, az, el]
for i = 0:t_test
    idx = find(traj(:,1) == i);
    if ~isempty(idx)
        cur_pt = traj(idx,:);
    end
    utc_p = datestr(datenum(utc, 'dd-mmm-yyyy HH:MM:SS') + minutes(i));
    [az, el] = radec2azel(cur_pt(3), cur_pt(2), lat, lon, utc_p);

    if el > el_min
        atraj_fine = [atraj_fine
                      i, az, el];
    else
        traj(idx,:) = [];
    end
end

fig = figure(fig_idx);

subplot(2,2,2);
hold on
plot(atraj_fine(:,1), rad2deg(atraj_fine(:,3)), 'b');
title('Elevation vs Time', 'FontSize', 9);
axis([0 t_test 0 90]);
line([0 t_test], rad2deg([el_min, el_min]), 'LineStyle', ':', 'Color', 'k');
xlabel('Time [min]', 'FontSize', 8);
ylabel('Elevation [deg]', 'FontSize', 8);
hold off


%% Histogram of elevations

fig = figure(fig_idx);

subplot(2,2,4);
histogram(rad2deg(atraj_fine(:,3)));    
title('Elevation Histogram', 'FontSize', 9);
ylabel('Time [min]', 'FontSize', 8);
xlabel('Elevation [deg]', 'FontSize', 8);


%% Final thoughts

fig = figure(fig_idx);

subplot(2, 2, 2);
mean_str = ['Mean: ' num2str(round(rad2deg(mean(atraj_fine(:,3))), 0)), '\circ'];
min_str = ['Min: ' num2str(round(rad2deg(min(atraj_fine(:,3))), 0)), '\circ'];

ax = gca;
offset = ax.XLim(2) * 0.02;
text(ax.XLim(2) + offset, 87, mean_str, 'FontSize', 8);
text(ax.XLim(2) + offset, 77, min_str, 'FontSize', 8);

% Outputs
% t = traj(:,1);  % Relative time
t = datenum(datestr(start_time + minutes(traj(:,1))));   % Absolute time
ra = rad2deg(traj(:,3));
dec = rad2deg(traj(:,2));

% Final cleanup
ra = wrapTo360(ra);

fprintf(['Done! (' num2str(toc) 's)' '\n'])