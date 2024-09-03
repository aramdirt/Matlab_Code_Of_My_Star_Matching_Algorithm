clear;clc;close all;

run ('my_config.m');
catalog_all = load('real_catalogue_with_sd.txt');
boresight = load('qtm_triangle_circumcenter.txt');

% Original point unit vector
FOV_O = zeros(1,3);
allstar_unitvector = zeros(length(catalog_all),3);
statistical_vector = zeros(1,100);
boresight_statistical_vector = zeros(length(boresight),4);

for i=1:1:length(boresight)
    
    RA_0 = boresight(i,7);
    DEC_0 = boresight(i,8);
    
    % Original FOV point unit vector
    FOV_O(1,1) = cosd(RA_0)*cosd(DEC_0);
    FOV_O(1,2) = sind(RA_0)*cosd(DEC_0);
    FOV_O(1,3) = sind(DEC_0);
    
    %all star unit vector
    for j=1:1:length(catalog_all)
        allstar_unitvector(j,1) = cosd(catalog_all(j,2))*cosd(catalog_all(j,3));
        allstar_unitvector(j,2)= sind(catalog_all(j,2))*cosd(catalog_all(j,3));
        allstar_unitvector(j,3) = sind(catalog_all(j,3));
    end
    
    %caculate number of star in FOV
     num_input_stars = 0;
    for j = 1:1:length(catalog_all)
        distance_all = rad2deg(acos(FOV_O(1,1)*allstar_unitvector(j,1) + FOV_O(1,2)*allstar_unitvector(j,2) + FOV_O(1,3)*allstar_unitvector(j,3)));
        if distance_all <= 6
            num_input_stars = num_input_stars + 1;
        end
    end
    statistical_vector(1,num_input_stars) = statistical_vector(1,num_input_stars) + 1;
    
    boresight_statistical_vector(i,1) = i;
    boresight_statistical_vector(i,2) = RA_0;
    boresight_statistical_vector(i,3) = DEC_0;
    boresight_statistical_vector(i,4) = num_input_stars;
end
boresight_statistical_vector = bubble_sort(boresight_statistical_vector,length(boresight_statistical_vector),4,1);
bar(statistical_vector);

% boresight_statistical_fd = fopen('boresight_statistical_data.txt','w+');
% fprintf(boresight_statistical_fd,'%d   %f   %f   %d \n', boresight_statistical_vector');
% fclose(boresight_statistical_fd);