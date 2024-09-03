clear;clc;close all;

run ('my_config.m');

%catalog_all = load('real_catalogue_without_doublestar.txt');
catalog_all = load('real_catalogue.txt');
%catalog_all = load('real_catalogue_without_doublestar_without_bright_star.txt');
%normal mode
catalog = load('guide_catalogue_1.txt');
catalog_ID = load('guide_catalogue_ID_1.txt');
mapping_table = load('mapping_table_1.txt'); 
%miss first nearby star mode
catalog_2 = load('guide_catalogue_2.txt');
catalog_ID_2 = load('guide_catalogue_ID_2.txt');
mapping_table_2 = load('mapping_table_2.txt');


RA_0=rand()*360;
DEC_0=rand()*180-90;
PHI_R=rand()*360;

RA_0= 35.048925;
DEC_0=  17.458732;
PHI_R=0;


limit_matching_number = 13;
limit_input_star_number = length(catalog_ID(1,:))-3;
% error_range = centroid_error*16; 
catalogue_FOV_size = 8.5;
Mv_error = 0;

% prepare input stars
F = FOCAL_LEN;
L = PIXEL_SIZE;
M = NUM_ROW; 
N = NUM_COL;
IMG_FOV_X = SYS_FOV_X;
IMG_FOV_Y = SYS_FOV_Y;
max_X = round((STR_FOV_X/IMG_FOV_X)*N/2);
max_Y = round((STR_FOV_Y/IMG_FOV_Y)*M/2);
frame_C_X = round(N/2);
X_end = frame_C_X+max_X;
frame_C_Y = round(M/2);
X_beg = frame_C_X-max_X;
Y_beg = frame_C_Y-max_Y;
Y_end = frame_C_Y+max_Y;

tic
% Original point unit vector
FOV_O = zeros(1,3);
[FOV_O(1,1),FOV_O(1,2),FOV_O(1,3)] = to_unit_vector(RA_0,DEC_0);

%all star unit vector
allstar_unitvector = zeros(length(catalog_all),3);
for i=1:1:length(catalog_all)
    [allstar_unitvector(i,1),allstar_unitvector(i,2),allstar_unitvector(i,3)] = to_unit_vector(catalog_all(i,2),catalog_all(i,3));
end

num_input_stars = 0;
for i = 1:length(catalog_all)
    distance_all = rad2deg(acos(FOV_O(1,1)*allstar_unitvector(i,1) + FOV_O(1,2)*allstar_unitvector(i,2) + FOV_O(1,3)*allstar_unitvector(i,3)));
    if distance_all <= 7.5
        num_input_stars = num_input_stars + 1;
        if num_input_stars > limit_input_star_number
            num_input_stars = limit_input_star_number;
            break;
        end
        input_stars(num_input_stars,1)=catalog_all(i,1);
        input_stars(num_input_stars,2)=catalog_all(i,2);
        input_stars(num_input_stars,3)=catalog_all(i,3);
        input_stars(num_input_stars,4)=catalog_all(i,4);    %bright
        input_stars(num_input_stars,5)=catalog_all(i,5);    %Smallest distance to nearby star(rad)
        input_stars(num_input_stars,6)=distance_all;         %distance to FOV center(deg)
    end
end



%choose main star
%input centroid error
for i=1:1:num_input_stars
    dec_ori=deg2rad(input_stars(i,3));
    
    error=rand()*centroid_error*2-centroid_error; % +- centroid_error
    
    ra_error=error/cos(dec_ori);  
    dec_error=error;
    bright_error = rand()*Mv_error*2-Mv_error;
    
    ra=deg2rad(input_stars(i,2)+ra_error/3600);
    dec=deg2rad(input_stars(i,3)+dec_error/3600);
    input_stars_unit_vector(i,1) = cos(ra)*cos(dec);
    input_stars_unit_vector(i,2) = sin(ra)*cos(dec);
    input_stars_unit_vector(i,3) = sin(dec);
    input_stars_unit_vector(i,4) = i; %which input star
    input_stars_unit_vector(i,5) = rad2deg(acos(FOV_O(1)*input_stars_unit_vector(i,1) + FOV_O(2)*input_stars_unit_vector(i,2) + FOV_O(3)*input_stars_unit_vector(i,3)));    %distance to FOV center(deg)
    input_stars_unit_vector(i,6) = input_stars(i,4)+bright_error;    %bright
end


input_stars_unit_vector = bubble_sort(input_stars_unit_vector,num_input_stars,5,1); %sort by (distance to FOV center)
mainstar(1) = input_stars_unit_vector(1,1); %mainstar unit vector X
mainstar(2) = input_stars_unit_vector(1,2); %mainstar unit vector Y
mainstar(3) = input_stars_unit_vector(1,3); %mainstar unit vector Z
mainstar(4) = input_stars(input_stars_unit_vector(1,4),1);%True mainstar ID
input_stars_unit_vector(1,:) = [];

%============use weight arrange nearby stars of mainstar================%
%PS : num_input_stars - 1 = number of nearby star
to_arrange_weight_nearby_star_number = 0;
for i=1:1:num_input_stars-1     % -1 because mainstar
    distance_to_mainstar = rad2deg(acos(mainstar(1)*input_stars_unit_vector(i,1) + mainstar(2)*input_stars_unit_vector(i,2) + mainstar(3)*input_stars_unit_vector(i,3))); %distance to mainstar(deg)
    %first delete far nearby star
    if distance_to_mainstar > catalogue_FOV_size
        continue;
    end
    to_arrange_weight_nearby_star_number = to_arrange_weight_nearby_star_number + 1;
    use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,1) = i;%which in (input_stars_unit_vector)
    use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,2) = input_stars_unit_vector(i,1);
    use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,3) = input_stars_unit_vector(i,2);
    use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,4) = input_stars_unit_vector(i,3);
    use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,5) = distance_to_mainstar;
    use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,6) = input_stars_unit_vector(i,6);%bright
end


[nearby_star_vector,input_count,ori_reliable_nearby_star_number] = use_weight_arrange_star(use_by_weight_arrange_vector,to_arrange_weight_nearby_star_number,input_stars,input_stars_unit_vector); 
% nearby_star_vector = [parterner_star;nearby_star_vector];

%PS : input_count here is number of reliable nearby star


%======================= LIS matching==========================%

[is_match_correct,num_matched_stars,final_output,num_mismatch,final_mainstar] = LIS_matching(input_count+1,mainstar,nearby_star_vector,catalog,catalog_ID,mapping_table,error_range,limit_matching_number); % +1 because mainstar

if is_match_correct == 0 &&  num_mismatch == 0
    [is_match_correct,num_matched_stars,final_output,num_mismatch,final_mainstar] = LIS_matching(input_count+1,mainstar,nearby_star_vector,catalog_2,catalog_ID_2,mapping_table_2,error_range,limit_matching_number);% +1 because mainstar
end
if is_match_correct == 0 &&  num_mismatch == 0
    nearby_star_vector(1,:) = [];
    input_count = input_count - 1;
    [is_match_correct,num_matched_stars,final_output,num_mismatch,final_mainstar] = LIS_matching(input_count+1,mainstar,nearby_star_vector,catalog,catalog_ID,mapping_table,error_range,limit_matching_number); % +1 because mainstar
end

z=inv(rand(1000));
toc

disp('Input:');    
disp(['   Original input stars numer = ',num2str(num_input_stars)]); 
disp(['   LIS input stars number = ',num2str(input_count + 1)]);    % +1 because mainstar
disp(['   Main Star = ', num2str(mainstar(4))]);
disp(['   Stars = ', num2str(nearby_star_vector(1:input_count,4)')]);

output_number = limit_matching_number;
if length(final_output) - 1 < limit_matching_number
    output_number = length(final_output) - 1;
end
disp('Final_Output');    
disp(['   Number of stars = ',num2str(num_matched_stars)]); 
disp(['   Main Star = ', num2str(final_mainstar)]);
disp(['   Match Stars = ', num2str(final_output(1,2:output_number+1))]);
disp(['   Number Mismatch Stars  = ', num2str(num_mismatch)]);

figure;
subplot(411);
plot(nearby_star_vector(1:input_count,4), '-*');
title(['Input Stars (RA=', num2str(RA_0), ', DEC=', num2str(DEC_0), ', PHI_R=', num2str(PHI_R), ')']);

subplot(412);
plot(final_output(1,2:input_count+1), '-ms');
title('final output stars ID');

for i = 1:1:input_count
    if final_output(1,i+1) == nearby_star_vector(i,4)
        results(i) = 1;
    else 
        results(i) = 0;
    end
end

subplot(413);
plot(results, '-ro');
title(['Match Results: ', num2str(num_matched_stars), '/', num2str(num_input_stars)]);
