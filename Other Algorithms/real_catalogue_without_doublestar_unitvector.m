clear;clc;close all;

catalog_all = load('real_catalogue_with_sd.txt');

for i=1:1:length(catalog_all)
    [catalog_all(i,6),catalog_all(i,7),catalog_all(i,8)] = to_unit_vector(catalog_all(i,2),catalog_all(i,3));
end
   
catalog_all_fd = fopen('real_catalogue.txt','w+');
fprintf(catalog_all_fd,'%d, %f, %f, %f, %f, %f, %f, %f\n', catalog_all');
fclose(catalog_all_fd);