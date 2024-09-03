clear;clc;close all;

run ('my_config.m');

catalog_ID = load('guide_catalogue_ID_1.txt');
catalog_ID_2 = load('guide_catalogue_ID_2.txt');

%TRK searching catalogue make 1
for i=1:1:length(catalog_ID(:,1))
    count = 0;
    for j=1:1:length(catalog_ID(1,:))
        if j == 2 || j == 3
            continue;
        end
        count = count + 1;
        for k =1:1:length(catalog_ID(:,1))
            if catalog_ID(i,j) == catalog_ID(k,1)
                if catalog_ID(i,j) == 0
                    catalog_TRK(i,count) = 0;
                else
                    catalog_TRK(i,count) = k;
                end
            end
        end
    end
end

%TRK searching catalogue make 2
for i=1:1:length(catalog_ID_2(:,1))
    count = 0;
    for j=1:1:length(catalog_ID_2(1,:))
        if j == 2 || j == 3
            continue;
        end
        count = count + 1;
        for k =1:1:length(catalog_ID_2(:,1))
            if catalog_ID_2(i,j) == catalog_ID_2(k,1)
                if catalog_ID_2(i,j) == 0
                    catalog_TRK_2(i,count) = 0;
                else
                    catalog_TRK_2(i,count) = k;
                end
            end
        end
    end
end

catalog_TRK_fd = fopen('catalogue_TRK_1.txt','w+');
for i=1:1:length(catalog_TRK(:,1))
    for j=1:1:length(catalog_TRK(1,:))
        fprintf(catalog_TRK_fd,'%d ', catalog_TRK(i,j));
    end
    fprintf(catalog_TRK_fd,'\n');
end
fclose(catalog_TRK_fd);

catalog_TRK_2_fd = fopen('catalogue_TRK_2.txt','w+');
for i=1:1:length(catalog_TRK_2(:,1))
    for j=1:1:length(catalog_TRK_2(1,:))
        fprintf(catalog_TRK_2_fd,'%d ', catalog_TRK_2(i,j));
    end
    fprintf(catalog_TRK_2_fd,'\n');
end
fclose(catalog_TRK_2_fd);

