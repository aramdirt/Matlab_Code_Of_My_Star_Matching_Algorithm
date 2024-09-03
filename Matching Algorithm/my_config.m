

%%%% Camera parameters
FOCAL_LEN = 16; % Focus length, UNIT: MM
PIXEL_SIZE = 2.2/1000; % Pixel size, UNIT: MM
NUM_ROW_PIXEL = 2592; % Image pixel number, number of pixel per row
NUM_COL_PIXEL = 1944; % Image pixel number, number of pixel per column
NUM_ROW = NUM_COL_PIXEL;    % number of rows
NUM_COL = NUM_ROW_PIXEL;    % number of columns
DEG_PER_PIXEL = 2*atan(PIXEL_SIZE/FOCAL_LEN/2)/pi*180; % DEG/PIXEL
SYS_FOV_X = NUM_ROW_PIXEL * DEG_PER_PIXEL; % X�����j�p
SYS_FOV_Y = NUM_COL_PIXEL * DEG_PER_PIXEL; % Y�����j�p
NOISE_LEVEL = 0;        % Noise level  
BACKGROUND_LEVEL = 50;  % Background level 

% FOV, degree
STR_FOV_X = 12;
STR_FOV_Y = 12;

% FOV Horizontal, degree
FOV_H=STR_FOV_X;

% FOV Diognal, degree
FOV_D=sqrt(STR_FOV_X^2+STR_FOV_Y^2);

% Maximum angular distance in degree
MAX_angular_dist=FOV_H;

% degree / pixel
deg_per_pixel=15.2249/1944;

% Minimum angular distance in degree
% MIN_angular_dist=5*deg_per_pixel;   % 5 pixels
MIN_angular_dist=rad2deg(0.0017);

% Unit for inter_star_distance, in arc-second
angular_dist_unit=(1/16);

% Unit for search-index table, in arc-second
idx_unit=1;

% raw key to idx key factor
mapping_factor=idx_unit/angular_dist_unit;

% mapping function for inter_star_distance to search index table
LIS_mapping_fun = @(ang_dist)(floor(ang_dist/mapping_factor)+1);

% error when centroiding, in arc-second
centroid_error=30;

% angular distance searching tolerance, in arc-second
if centroid_error~=0
    error_range=(centroid_error*2+5)*16;
%       searching_tolerance=12;
else
    error_range=1*16;
end

% Maximum matching stars
Max_matching_stars=15;

% Minimum required stars for matching
Min_matching_stars=4;

% STR attitude output rate, in Hz
STR_output_rate=4;

% STR maximum moving per axis, degree/second
Max_STR_moving_per_axis=2;

% use small or large block_inter_star_distance?
use_large_block_inter_star_distance=0;

% use star's Mv in matching?
use_Mv_in_gva_voting=0;

