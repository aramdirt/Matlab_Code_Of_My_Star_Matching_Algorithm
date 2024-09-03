clear;clc;close all;

Catalog_1 = load('Guide_Catalog_1.txt');
Catalog_2 = load('Guide_Catalog_2.txt');
Searching_List_2 = load('Searching_List_2.txt');

for i=1:1:length(Catalog_1(:,1))
    Guide_Catalogue(i,:) = Catalog_1(i,:);
%     for j=1:1:length(Catalog_2(:,1))
%         if Catalog_1(i,1) == Catalog_2(j,1)
%             Guide_Catalogue(i,9) = Catalog_2(j,8);
%         end
%     end
end

guide_catalog_fd = fopen('Guide_Catalogue.txt','w+');
for i=1:1:length(Guide_Catalogue(:,1))
    for j=1:1:length(Guide_Catalogue(1,:))
        if j == length(Guide_Catalogue(1,:))
            fprintf(guide_catalog_fd,'%s', num2str(Guide_Catalogue(i,j) - 1));
            fprintf(guide_catalog_fd,'\n');
            break;
        end
        fprintf(guide_catalog_fd,'%s, ', num2str(Guide_Catalogue(i,j)));
    end
end
fclose(guide_catalog_fd);