clear all;clc;close all;

real_catalog = load('real_catalogue_without_doublestar.txt');
count = 0;
for i=1:1:length(real_catalog)
    if real_catalog(i,4) >= 1
        count = count + 1;
        catalog(count,:) = real_catalog(i,:);
    end
end