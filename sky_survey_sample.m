% Sample script. Generates an arbitrary sky survey

clear
clc

% Physical parameters
el_min = 40;
t_test = 5 * 60;
fov = 7.5;
lat = 23;
lon = -79;

% Tuning parameters
tuning = struct('prioritize_high_elevations', true, ...
            'timing', ...
                struct('adjust', true, ...
                       'min_step', 1), ...
            'proximity', ...
                struct('adjust', true, ...
                      'window_cost_threshold', 0.3, ...
                      'window_max_size', 10, ...
                      'ra_weight', 1, ...
                      'dec_weight', 1));

% Misc parameters
fig_idx = 50;
show_elevs = true;

% Go time
[t, ra, dec, fig] = sky_survey(el_min, t_test, fov, lat, lon, tuning, fig_idx, show_elevs);