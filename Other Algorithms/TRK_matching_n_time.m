clear;clc;close all;

run ('my_config.m');

catalog_all = load('real_catalogue_with_sd.txt');
%normal mode
catalog = load('guide_catalogue_1.txt');
catalog_ID = load('guide_catalogue_ID_1.txt');
catalog_TRK = load('catalogue_TRK_1.txt');
%miss first nearby star mode
catalog_2 = load('guide_catalogue_2.txt');
catalog_ID_2 = load('guide_catalogue_ID_2.txt');
catalog_TRK_2 = load('catalogue_TRK_2.txt');

limit_matching_number = 15;
limit_input_star_number = 109;
error_range = 0.01;

loop_count = 100;
total_match_star = 0;
total_mismatch_star = 0;
successful_test = 0;

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
frame_C_Y = round(M/2);
X_beg = frame_C_X-max_X;
X_end = frame_C_X+max_X;
Y_beg = frame_C_Y-max_Y;
Y_end = frame_C_Y+max_Y;

%all star unit vector
allstar_unitvector = zeros(length(catalog_all),3);
for i=1:1:length(catalog_all)
    [allstar_unitvector(i,1),allstar_unitvector(i,2),allstar_unitvector(i,3)] = to_unit_vector(catalog_all(i,2),catalog_all(i,3));
end
%loop start
for k = 1:1:loop_count
    clear dis2FOV;
    clear initial_nearbystar;
    clear searching_catalogue;
    clear searching_catalogue_ID;
    clear input_stars;
    clear input_stars_unit_vector;
    clear use_by_weight_arrange_vector;
    clear final_output;
    
    
    % create initial FOV
    RA_0=rand()*360;
    DEC_0=rand()*180-90;
    PHI_R=rand()*360;

    I_FOV_O = zeros(1,3);
    [I_FOV_O(1,1),I_FOV_O(1,2),I_FOV_O(1,3)] = to_unit_vector(RA_0,DEC_0);

    % create initial mainstar
    dis2FOV = zeros(length(catalog_all),1);
    for i =1:1:length(catalog_all)
        dis2FOV(i,1) = i;
        dis2FOV(i,2) =  acos(I_FOV_O(1,1)*allstar_unitvector(i,1) + I_FOV_O(1,2)*allstar_unitvector(i,2) + I_FOV_O(1,3)*allstar_unitvector(i,3));
    end
    dis2FOV = bubble_sort(dis2FOV,length(dis2FOV),2,1);

    for i=1:1:length(catalog_ID(:,1))
        if dis2FOV(1,1) == catalog_ID(i,1)
            initial_mainstar = i; %position in catalog_ID
            break;
        end
    end
    initial_nearbystar = catalog_TRK(initial_mainstar,2:length(catalog_TRK(1,:)));  %position in catalog_ID



    % STR Moving
    tracking_moving_dist=Max_STR_moving_per_axis/STR_output_rate;
    ra_moving=rand()*tracking_moving_dist*2-tracking_moving_dist;
    ra_moving=ra_moving/cos(deg2rad(DEC_0));
    dec_moving=rand()*tracking_moving_dist*2-tracking_moving_dist;
    phi_moving=rand()*tracking_moving_dist*2-tracking_moving_dist;

    RA_0=abs(rem(RA_0+ra_moving,360));
    DEC_0=DEC_0+dec_moving;
    if DEC_0 >=90
        DEC_0=90-rem(DEC_0,90);
    end
    if DEC_0 <= -90
        DEC_0=-90-rem(DEC_0,90);
    end
    PHI_R=abs(rem(PHI_R+phi_moving,360));

    %make searching catalogue
    searching_catalogue(1,:) = catalog(initial_mainstar,:);
    searching_catalogue_ID(1,:) = catalog_ID(initial_mainstar,:);
    count = 1;
    for i=1:1:length(initial_nearbystar)
        if initial_nearbystar(i) == 0
            continue;
        end
        count = count + 1;
        searching_catalogue(count,:) = catalog(initial_nearbystar(i),:);
        searching_catalogue_ID(count,:) = catalog_ID(initial_nearbystar(i),:);
    end

    % NEW FOV center unit vector
    FOV_O = zeros(1,3);
    [FOV_O(1,1),FOV_O(1,2),FOV_O(1,3)] = to_unit_vector(RA_0,DEC_0);

    % input star to NEW  FOV
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
    
        ra=deg2rad(input_stars(i,2)+ra_error/3600);
        dec=deg2rad(input_stars(i,3)+dec_error/3600);
        input_stars_unit_vector(i,1) = cos(ra)*cos(dec);
        input_stars_unit_vector(i,2) = sin(ra)*cos(dec);
        input_stars_unit_vector(i,3) = sin(dec);
        input_stars_unit_vector(i,4) = i; %which input star
        input_stars_unit_vector(i,5) = rad2deg(acos(FOV_O(1)*input_stars_unit_vector(i,1) + FOV_O(2)*input_stars_unit_vector(i,2) + FOV_O(3)*input_stars_unit_vector(i,3)));    %distance to FOV center(deg)
        input_stars_unit_vector(i,6) = input_stars(i,4);    %bright
    end

    % define new mainstar
    input_stars_unit_vector = bubble_sort(input_stars_unit_vector,num_input_stars,5,1); %sort by (distance to FOV center)
    mainstar(1) = input_stars_unit_vector(1,1); %mainstar unit vector X
    mainstar(2) = input_stars_unit_vector(1,2); %mainstar unit vector Y
    mainstar(3) = input_stars_unit_vector(1,3); %mainstar unit vector Z
    mainstar(4) = input_stars(input_stars_unit_vector(1,4),1);%True mainstar ID


    %============use weight arrange nearby stars of mainstar================%
    %PS : num_input_stars - 1 = number of nearby star
    to_arrange_weight_nearby_star_number = 0;
    for i=1:1:num_input_stars-1     % -1 because mainstar
        distance_to_mainstar = rad2deg(acos(mainstar(1)*input_stars_unit_vector(i+1,1) + mainstar(2)*input_stars_unit_vector(i+1,2) + mainstar(3)*input_stars_unit_vector(i+1,3))); %distance to mainstar(deg)
        %first delete far nearby star
        if distance_to_mainstar > 8.5
            continue;
        end
        to_arrange_weight_nearby_star_number = to_arrange_weight_nearby_star_number + 1;
        use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,1) = i+1;%which in (input_stars_unit_vector)
        use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,2) = input_stars_unit_vector(i+1,1);
        use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,3) = input_stars_unit_vector(i+1,2);
        use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,4) = input_stars_unit_vector(i+1,3);
        use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,5) = distance_to_mainstar;
        use_by_weight_arrange_vector(to_arrange_weight_nearby_star_number,6) = input_stars_unit_vector(i+1,6);%bright
    end

    [nearby_star_vector,input_count,ori_reliable_nearby_star_number] = use_weight_arrange_star(use_by_weight_arrange_vector,to_arrange_weight_nearby_star_number,input_stars,input_stars_unit_vector); 
    %PS : input_count here is number of reliable nearby star

    %======================= TRK matching==========================%
    [is_match_correct,num_matched_stars,final_output,num_mismatch,final_mainstar] = TRK_matching(input_count+1,mainstar,nearby_star_vector,searching_catalogue,searching_catalogue_ID,error_range,limit_matching_number); % +1 because mainstar
    
    
    total_match_star = total_match_star + num_matched_stars;
    total_mismatch_star = total_mismatch_star + num_mismatch;
    if is_match_correct == 1
        successful_test = successful_test + 1;
    end

end
disp(['   total match star = ',num2str(total_match_star)]); 
disp(['   average match star = ',num2str(total_match_star/loop_count)]); 
disp(['   total mismatch star = ',num2str(total_mismatch_star)]);
disp(['   successful test = ',num2str(successful_test)]);
