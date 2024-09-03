clear;clc;close all;

catalog = load('guide_catalogue_1.txt');
catalog_ID = load('guide_catalogue_ID_1.txt');
catalog_TRK = load('catalogue_TRK_1.txt');
catalog_mapping = load('mapping_table_1.txt');

catalog_2 = load('guide_catalogue_2.txt');
catalog_ID_2 = load('guide_catalogue_ID_2.txt');
catalog_TRK_2 = load('catalogue_TRK_2.txt');
catalog_mapping_2 = load('mapping_table_2.txt');

real_catalog = load('real_catalogue.txt');

%% Catalogue 1 %%
guide_catalog_1 = zeros(length(catalog(:,1)),7);
searching_list_count_1 = 0;

for i=1:1:length(catalog(:,1))
    %guide_catalog_1(i,1) = i;   %idx
    guide_catalog_1(i,1) = catalog(i,1); %Main Star ID
    guide_catalog_1(i,2) = real_catalog(catalog(i,1),2);   %RA
    guide_catalog_1(i,3) = real_catalog(catalog(i,1),3);    %DEC
    guide_catalog_1(i,4) = catalog(i,2); %distance to Partner Star (angular distance unit)
    guide_catalog_1(i,5) = catalog_ID(i,4); %Partner Star ID
    guide_catalog_1(i,6) = catalog_TRK(i,2) - 1; %Partner Star idx
    guide_catalog_1(i,7) = catalog(i,3)-1; %Nearby Star Number
    if i ~= 1
        guide_catalog_1(i,8) = guide_catalog_1(i-1,7) + guide_catalog_1(i-1,8);   %searching head
    else
        guide_catalog_1(i,8) = 1;   %searching head
    end
    searching_list_count_1 = searching_list_count_1 + guide_catalog_1(i,7);
    count = 1;
    for j=guide_catalog_1(i,8):1:(guide_catalog_1(i,8) + guide_catalog_1(i,7) - 1)
        searching_list_1(j,1) = catalog_ID(i,count+4);  %nearby star ID
        searching_list_1(j,2) = catalog(i,(2*count)+2);    %feature 1
        searching_list_1(j,3) = catalog(i,(2*count)+3);    %feature 2
        searching_list_1(j,4) = catalog_TRK(i,count+2) - 1;    %nearby star idx 
        count = count + 1;
    end
end 

guide_catalog_1_fd = fopen('Guide_Catalog_1.txt','w+');
for i=1:1:length(guide_catalog_1(:,1))
    for j=1:1:length(guide_catalog_1(1,:))
        if j == length(guide_catalog_1(1,:))
            fprintf(guide_catalog_1_fd,'%s', num2str(guide_catalog_1(i,j) - 1));
            fprintf(guide_catalog_1_fd,'\n');
            break;
        end
        fprintf(guide_catalog_1_fd,'%s, ', num2str(guide_catalog_1(i,j)));
    end
end
fclose(guide_catalog_1_fd);

searching_list_1_fd = fopen('Searching_List_1.txt','w+');
for i=1:1:length(searching_list_1(:,1))
    for j=1:1:length(searching_list_1(1,:))
        if j == length(searching_list_1(1,:))
            fprintf(searching_list_1_fd,'%s', num2str(searching_list_1(i,j)));
            fprintf(searching_list_1_fd,'\n');
            break;
        end
        fprintf(searching_list_1_fd,'%s, ', num2str(searching_list_1(i,j)));
    end
end
fclose(searching_list_1_fd);

%% Mapping Table 1 %%
mapping_talbe_1_fd = fopen('Mapping_Table_1_C.txt','w+');
for i=1:1:length(catalog_mapping(:,1))
    fprintf(mapping_talbe_1_fd,'%s,%s\n', num2str(catalog_mapping(i,2)),num2str(catalog_mapping(i,3)));
end
fclose(mapping_talbe_1_fd);

%% Catalogue 2 %%

guide_catalog_2 = zeros(length(catalog_2(:,1)),7);
searching_list_count_2 = 0;

for i=1:1:length(catalog_2(:,1))
    %guide_catalog_2(i,1) = i;   %idx
    guide_catalog_2(i,1) = catalog_2(i,1); %Main Star ID
    guide_catalog_2(i,2) = real_catalog(catalog_2(i,1),2);   %RA
    guide_catalog_2(i,3) = real_catalog(catalog_2(i,1),3);    %DEC
    guide_catalog_2(i,4) = catalog_2(i,2); %distance to Partner Star (angular distance unit)
    guide_catalog_2(i,5) = catalog_ID_2(i,4); %Partner Star ID
    guide_catalog_2(i,6) = catalog_TRK_2(i,2) - 1; %Partner Star idx
    guide_catalog_2(i,7) = catalog_2(i,3)-1; %Nearby Star Number
    if i ~= 1
        guide_catalog_2(i,8) = guide_catalog_2(i-1,7) + guide_catalog_2(i-1,8);   %searching head
    else
        guide_catalog_2(i,8) = 1;   %searching head
    end
    searching_list_count_2 = searching_list_count_2 + guide_catalog_2(i,7);
    count = 1;
    for j=guide_catalog_2(i,8):1:(guide_catalog_2(i,8) + guide_catalog_2(i,7) - 1)
        searching_list_2(j,1) = catalog_ID_2(i,count+4);  %nearby star ID
        searching_list_2(j,2) = catalog_2(i,(2*count)+2);    %feature 1
        searching_list_2(j,3) = catalog_2(i,(2*count)+3);    %feature 2
        searching_list_2(j,4) = catalog_TRK_2(i,count+2) - 1;    %nearby star idx
        count = count + 1;
    end
end 

guide_catalog_2_fd = fopen('Guide_Catalog_2.txt','w+');
for i=1:1:length(guide_catalog_2(:,1))
    for j=1:1:length(guide_catalog_2(1,:))
        if j == length(guide_catalog_2(1,:))
            fprintf(guide_catalog_2_fd,'%s', num2str(guide_catalog_2(i,j) - 1));
            fprintf(guide_catalog_2_fd,'\n');
            break;
        end
        fprintf(guide_catalog_2_fd,'%s, ', num2str(guide_catalog_2(i,j)));
    end
end
fclose(guide_catalog_2_fd);

searching_list_2_fd = fopen('Searching_List_2.txt','w+');
for i=1:1:length(searching_list_2(:,1))
    for j=1:1:length(searching_list_2(1,:))
        if j == length(searching_list_2(1,:))
            fprintf(searching_list_2_fd,'%s', num2str(searching_list_2(i,j)));
            fprintf(searching_list_2_fd,'\n');
            break;
        end
        fprintf(searching_list_2_fd,'%s, ', num2str(searching_list_2(i,j)));
    end
end
fclose(searching_list_2_fd);

%% Mapping Table 2 %%
mapping_talbe_2_fd = fopen('Mapping_Table_2_C.txt','w+');
for i=1:1:length(catalog_mapping_2(:,1))
    fprintf(mapping_talbe_2_fd,'%s,%s\n', num2str(catalog_mapping_2(i,2)),num2str(catalog_mapping_2(i,3)));
end
fclose(mapping_talbe_2_fd);