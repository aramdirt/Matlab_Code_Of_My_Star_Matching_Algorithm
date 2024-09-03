clear;clc;close all;

run ('my_config.m');

catalog_all = load('real_catalogue_without_doublestar.txt');
%normal mode
catalog = load('guide_catalogue_1.txt');
catalog_ID = load('guide_catalogue_ID_1.txt');
mapping_table = load('mapping_table_1.txt'); 
%miss first nearby star mode
catalog_2 = load('guide_catalogue_2.txt');
catalog_ID_2 = load('guide_catalogue_ID_2.txt');
mapping_table_2 = load('mapping_table_2.txt');

boresight = load('qtm_triangle_circumcenter.txt');

RA_0 = 0;
DEC_0 = 0;
PHI_R = 0;

limit_matching_number = 13;
limit_input_star_number = length(catalog_ID(1,:))-3;
% error_range = searching_tolerance*16;
catalogue_FOV_size = 8.5;
Mv_error = 0;

F = FOCAL_LEN;
L = PIXEL_SIZE;
M = NUM_ROW; 
N = NUM_COL;
IMG_FOV_X = SYS_FOV_X;
IMG_FOV_Y = SYS_FOV_Y;
max_X = round((STR_FOV_X/IMG_FOV_X)*N/2);
max_Y = round((STR_FOV_Y/IMG_FOV_Y)*M/2);
frame_C_X = round(N/2);
frame_C_Y = round(M/2);
X_beg = frame_C_X-max_X;
X_end = frame_C_X+max_X;
Y_beg = frame_C_Y-max_Y;
Y_end = frame_C_Y+max_Y;    

% Original point unit vector
FOV_O = zeros(1,3);
allstar_unitvector = zeros(length(catalog_all),3);

%all star unit vector
allstar_unitvector = zeros(length(catalog_all),3);
for i=1:1:length(catalog_all)
    [allstar_unitvector(i,1),allstar_unitvector(i,2),allstar_unitvector(i,3)] = to_unit_vector(catalog_all(i,2),catalog_all(i,3));
end

%================Loop start=================%
num_tests = 0;
NG_cnts=0;
num_incorrect=0;

for for_counter = 1:1:length(boresight)
    clear input_stars;
    clear input_stars_unit_vector;
    clear nearby_star_vector;
    clear use_by_weight_arrange_vector;
    clear final_output;
    
    RA_0 = boresight(for_counter,7);
    DEC_0 = boresight(for_counter,8);
    
    % Original point unit vector
    FOV_O = zeros(1,3);
    [FOV_O(1,1),FOV_O(1,2),FOV_O(1,3)] = to_unit_vector(RA_0,DEC_0);
    

    
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
            input_stars(num_input_stars,6)=distance_all;         %distance to FOV center
        end
    end

    
    %choose main star
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
    input_stars_unit_vector = bubble_sort(input_stars_unit_vector,num_input_stars,5,1);
    mainstar(1) = input_stars_unit_vector(1,1); %mainstar unit vector X
    mainstar(2) = input_stars_unit_vector(1,2); %mainstar unit vector Y
    mainstar(3) = input_stars_unit_vector(1,3); %mainstar unit vector Z
    mainstar(4) = input_stars(input_stars_unit_vector(1,4),1);%True mainstar ID
    
    input_stars_unit_vector(1,:) = [];

    
    %use weight arrange other stars
    %PS : num_input_stars - 1 = number of nearby star
    to_arrange_weight_nearby_star_number = 0;
    for i=1:1:num_input_stars-1 % -1 because mainstar
        distance_to_mainstar = rad2deg(acos(mainstar(1)*input_stars_unit_vector(i,1) + mainstar(2)*input_stars_unit_vector(i,2) + mainstar(3)*input_stars_unit_vector(i,3))); %distance to mainstar(deg)
        if distance_to_mainstar > catalogue_FOV_size
            continue;
        end
        to_arrange_weight_nearby_star_number = to_arrange_weight_nearby_star_number + 1;
        use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,1) = i;%which in input_stars_unit_vector
        use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,2) = input_stars_unit_vector(i,1);
        use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,3) = input_stars_unit_vector(i,2);
        use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,4) = input_stars_unit_vector(i,3);
        use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,5) = distance_to_mainstar;
        use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,6) = input_stars_unit_vector(i,6);%bright
    end
    
    [nearby_star_vector,input_count,ori_reliable_nearby_star_number] = use_weight_arrange_star(use_by_weight_arrange_vector,to_arrange_weight_nearby_star_number,input_stars,input_stars_unit_vector);
%     nearby_star_vector = [parterner_star;nearby_star_vector];

    %======================= LIS matching==========================%
    [is_match_correct,num_matched_stars,final_output,num_mismatch,final_mainstar] = LIS_matching(input_count+1,mainstar,nearby_star_vector,catalog,catalog_ID,mapping_table,error_range,limit_matching_number);
    if is_match_correct == 0 &&  num_mismatch == 0
        [is_match_correct,num_matched_stars,final_output,num_mismatch,final_mainstar] = LIS_matching(input_count+1,mainstar,nearby_star_vector,catalog_2,catalog_ID_2,mapping_table_2,error_range,limit_matching_number);
    end
    if is_match_correct == 0 &&  num_mismatch == 0
        nearby_star_vector(1,:) = [];
        input_count = input_count - 1;
        [is_match_correct,num_matched_stars,final_output,num_mismatch,final_mainstar] = LIS_matching(input_count+1,mainstar,nearby_star_vector,catalog,catalog_ID,mapping_table,error_range,limit_matching_number); % +1 because mainstar
    end
    
    num_tests = num_tests+1;
    test_result(num_tests,1)=num_tests;
    test_result(num_tests,2)=num_matched_stars;
    test_result(num_tests,3)=RA_0;
    test_result(num_tests,4)=DEC_0;
    test_result(num_tests,5)=PHI_R;
    
    if is_match_correct ~= 1
        num_incorrect=num_incorrect+1;
    end
    
    if num_matched_stars +  1 < Min_matching_stars || ~is_match_correct     %num_matched_stars + 1 is mainstar
        NG_cnts=NG_cnts+1;
        NG_result(NG_cnts,1)=num_tests;
        if is_match_correct == 1
            NG_result(NG_cnts,2)=num_matched_stars + 1;
        else
            NG_result(NG_cnts,2)=num_matched_stars;
        end
        NG_result(NG_cnts,3)=num_mismatch;
        NG_result(NG_cnts,4)=RA_0;
        NG_result(NG_cnts,5)=DEC_0;
        NG_result(NG_cnts,6)=PHI_R;
    end
    disp(['   for_counter = ',num2str(for_counter)]);

end

if NG_cnts==0
    NG_result=zeros(1,5);
end
 NG_fd = fopen('LIS_coverage_NG_30.txt','w+');
fprintf(NG_fd,'%4d, %2d, %2d, %.8f, %.8f, %.8f\n',NG_result');
fclose(NG_fd);

disp(['Total tests = ',num2str(num_tests)]);
disp('Matching Correctness:');
disp(['   Incorrect tests = ',num2str(num_incorrect)]);
disp(['   Correct Ratio = ',num2str((1-num_incorrect/num_tests)*100),' %']);

disp('Matching Performance:');
disp(['   Correct tests = ',num2str(num_tests-num_incorrect)]);
disp(['   Matched stars < 4 tests = ',num2str(NG_cnts-num_incorrect)]);
disp(['   Successful Ratio = ',num2str((1-(NG_cnts-num_incorrect)/(num_tests-num_incorrect))*100),' %']);