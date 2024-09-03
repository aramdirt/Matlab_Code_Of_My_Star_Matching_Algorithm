clear;clc;close all;
catalog_all = load('real_catalogue_catST16_J2017.txt');

star_unit_vector = zeros(length(catalog_all),3);
for i=1:1:length(catalog_all)
    star_unit_vector(i,1) = cosd(catalog_all(i,2))*cosd(catalog_all(i,3));
    star_unit_vector(i,2)= sind(catalog_all(i,2))*cosd(catalog_all(i,3));
    star_unit_vector(i,3) = sind(catalog_all(i,3));
end
smallset_dis = zeros(length(catalog_all),1);
for i=1:1:length(catalog_all)
    clear distance;
    for j=1:1:length(catalog_all)
        distance(j,1) = acos(star_unit_vector(i,1)*star_unit_vector(j,1) + star_unit_vector(i,2)*star_unit_vector(j,2) + star_unit_vector(i,3)*star_unit_vector(j,3));
        if i == j
            distance(j,1) = 100;
        end
    end
    distance = bubble_sort(distance,length(catalog_all),1,1);
    catalog_all(i,5) = distance(1,1);
end

