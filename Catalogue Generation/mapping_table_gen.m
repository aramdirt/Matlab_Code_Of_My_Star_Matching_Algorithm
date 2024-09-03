clear;clc;close all;

run ('my_config.m');
catalog = load('guide_catalogue_1.txt');
catalog_2 = load('guide_catalogue_2.txt');
Max_Catalogue_star_number = length(catalog(:,1));
one_rad_angular_distance_unit = 33005;  %angular distance unit(0.01 rad)
constrain = one_rad_angular_distance_unit;  
count = 0;
for i=1:1:length(catalog(:,1))
    if catalog(i,2) >= constrain
        constrain = constrain + one_rad_angular_distance_unit;
        count = count + 1;
        sort(1,count) = i;
    end
end


for i=1:1:length(sort)
    mapping_table(i,1) = i;
    if i - 1 <= 0
        mapping_table(i,2) = 1;
    else
        mapping_table(i,2) = sort(1,i-1);
    end
    
    if i + 1 > floor(catalog(length(catalog(:,1)),2)/one_rad_angular_distance_unit)
        mapping_table(i,3) = Max_Catalogue_star_number;
    else
        mapping_table(i,3) = sort(1,i+1)-1;
    end
end

mapping_table_fd = fopen('mapping_table_1.txt','w+');
fprintf(mapping_table_fd,'%d ,%d ,%d \n', mapping_table');
fclose(mapping_table_fd);
constrain = one_rad_angular_distance_unit;  %angular distance unit(0.01 rad)
count_2 = 0;
%second mapping table
for i=1:1:length(catalog_2(:,1))
    if catalog_2(i,2) >= constrain
        constrain = constrain + one_rad_angular_distance_unit;
        count_2 = count_2 + 1;
        sort_2(1,count_2) = i;
    end
end


for i=1:1:length(sort_2)
    mapping_table_2(i,1) = i;
    if i - 1 <= 0
        mapping_table_2(i,2) = 1;
    else
        mapping_table_2(i,2) = sort_2(1,i-1);
    end
    
    if i + 1 > floor(catalog_2(length(catalog_2(:,1)),2)/one_rad_angular_distance_unit)
        mapping_table_2(i,3) = Max_Catalogue_star_number;
    else
        mapping_table_2(i,3) = sort_2(1,i+1)-1;
    end
end

mapping_table_2_fd = fopen('mapping_table_2.txt','w+');
fprintf(mapping_table_2_fd,'%d ,%d ,%d \n', mapping_table_2');
fclose(mapping_table_2_fd);