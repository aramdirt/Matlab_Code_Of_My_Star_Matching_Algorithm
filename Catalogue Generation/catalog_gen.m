clear all;clc;close all;
run ('my_config.m');

catalog = load('real_catalogue_without_doublestar.txt');
RA_idx = 2;
DEC_idx = 3;
smallest_input_number = 4;
num_stars_in_FOV=zeros(length(catalog),1);
FOV_radius = 8.5;


for i=1:1:length(catalog)
    clear stars_in_FOV_ID; 
    clear stars_in_FOV_CSS;
    clear fail_CSS_star;
    clear OK_CSS_star;
    clear stars_in_FOV;
    fail_count = 0;
    OK_count = 0;
    
    RA_0 = catalog(i,RA_idx);
    DEC_0 = catalog(i,DEC_idx);
    PHI_R = 0;
    [Main_star_vector(1),Main_star_vector(2),Main_star_vector(3)] = to_unit_vector(catalog(i,RA_idx),catalog(i,DEC_idx));
    
    [star_x_i,star_y_i,star_z_i] = cel2xy(catalog(i,RA_idx),catalog(i,DEC_idx),RA_0,DEC_0,PHI_R);

%%  for debug
    if i==5424
        C = 0;
    end
%%
    
    for j=1:1:length(catalog)
        if i==j
            continue;
        end
        
        [star_x_j,star_y_j,star_z_j] = cel2xy(catalog(j,RA_idx),catalog(j,DEC_idx),RA_0,DEC_0,PHI_R);
      
        angular_dist=acos(star_x_i*star_x_j+star_y_i*star_y_j+star_z_i*star_z_j); % in radian
        angular_dist=rad2deg(angular_dist); % in degree
        
    if angular_dist <= FOV_radius
        num_stars_in_FOV(i) = num_stars_in_FOV(i)+1;
        idx = num_stars_in_FOV(i);
        stars_in_FOV(idx,1) = j;
        stars_in_FOV(idx,2) = catalog(j,5);     %smallest destance to nearby star(rad)
        stars_in_FOV(idx,3) = angular_dist;     %distance to mainstar(deg)
        stars_in_FOV(idx,4) = catalog(j,4);     %star bright
    end
    end
    
    
    for k=1:1:length(stars_in_FOV(:,1))
        stars_in_FOV_CSS(k,1) = stars_in_FOV(k,1);
        [stars_in_FOV_CSS(k,2),stars_in_FOV_CSS(k,3)] = to_weight(stars_in_FOV(k,2),stars_in_FOV(k,3));
        stars_in_FOV_CSS(k,4) = stars_in_FOV(k,3);  %distance to mainstar(deg)
        stars_in_FOV_CSS(k,5) = stars_in_FOV(k,4);  %star bright
    end
    
    
    
    %CSS weight of all input star
    stars_in_FOV_CSS = bubble_sort(stars_in_FOV_CSS,length(stars_in_FOV_CSS(:,1)),3,1);
    for k=1:1:length(stars_in_FOV_CSS(:,1))
        if stars_in_FOV_CSS(k,3) == 0       %this star is double star with mainstar
            continue;
        end
        if stars_in_FOV_CSS(k,2) == 0
            fail_count = fail_count + 1;
            fail_CSS_star(fail_count,:) = stars_in_FOV_CSS(k,:);
        else
            OK_count = OK_count + 1;
            OK_CSS_star(OK_count,:)  = stars_in_FOV_CSS(k,:);
        end
    end
    %deal with OK CSS star first
    OK_CSS_star = bubble_sort(OK_CSS_star,OK_count,3,0);    %big to small
    first_flag = 0;
    input_count = 0;        % final input star number

    
    if OK_count > 1
        for k=2:1:OK_count
            if first_flag == 0
               input_count = input_count + 1;
               stars_in_FOV_ID(input_count,1) = OK_CSS_star(1,1);
               stars_in_FOV_ID(input_count,2) = OK_CSS_star(1,4);
               first_flag = 1;
            end
            % if same weight, arranged by bright
            if OK_CSS_star(k,3) + 0.008 >= OK_CSS_star(k-1,3)
                if OK_CSS_star(k,5) < OK_CSS_star(k-1,5)
                    input_count = input_count + 1;
                    stars_in_FOV_ID(input_count,1) = stars_in_FOV_ID(input_count - 1,1);
                    stars_in_FOV_ID(input_count,2) = stars_in_FOV_ID(input_count - 1,2);
                    
                    stars_in_FOV_ID(input_count - 1,1) = OK_CSS_star(k,1);      
                    stars_in_FOV_ID(input_count - 1,2) = OK_CSS_star(k,4);
                else
                    input_count = input_count + 1;
                    stars_in_FOV_ID(input_count,1) = OK_CSS_star(k,1);
                    stars_in_FOV_ID(input_count,2) = OK_CSS_star(k,4);
                end
            else
            input_count = input_count + 1;
            stars_in_FOV_ID(input_count,1) = OK_CSS_star(k,1);
            stars_in_FOV_ID(input_count,2) = OK_CSS_star(k,4);
            end
        end
    else
        input_count = input_count + 1;
        stars_in_FOV_ID(input_count,1) = OK_CSS_star(1,1);
        stars_in_FOV_ID(input_count,2) = OK_CSS_star(1,4);
    end

    % to unit vector
    for m=1:1:input_count
        stars_in_FOV_ID(m,3) = cosd(catalog(stars_in_FOV_ID(m,1),RA_idx))*cosd(catalog(stars_in_FOV_ID(m,1),DEC_idx));%X
        stars_in_FOV_ID(m,4) = sind(catalog(stars_in_FOV_ID(m,1),RA_idx))*cosd(catalog(stars_in_FOV_ID(m,1),DEC_idx));%Y
        stars_in_FOV_ID(m,5) = sind(catalog(stars_in_FOV_ID(m,1),DEC_idx));%Z
    end
    %=========================normal mode============================%
    %SPD
    %normal mode
    for k=1:1:(input_count-1)
        VB(i,1) = catalog(i,1);
        VB(i,2) = (rad2deg(acos(stars_in_FOV_ID(1,3)*Main_star_vector(1) + stars_in_FOV_ID(1,4)*Main_star_vector(2) + stars_in_FOV_ID(1,5)*Main_star_vector(3)))*3600)*16; % distance to first (angular distance unit)
        VB(i,3) = input_count;
        VB(i,(2*k) + 2) = (rad2deg(acos(stars_in_FOV_ID((k+1),3)*stars_in_FOV_ID(1,3) + stars_in_FOV_ID((k+1),4)*stars_in_FOV_ID(1,4) + stars_in_FOV_ID((k+1),5)*stars_in_FOV_ID(1,5)))*3600)*16;  %nearby 2 partner (angular distance unit)
        VB(i,(2*k) + 3) = (rad2deg(acos(stars_in_FOV_ID((k+1),3)* Main_star_vector(1) + stars_in_FOV_ID((k+1),4)* Main_star_vector(2) + stars_in_FOV_ID((k+1),5)* Main_star_vector(3)))*3600)*16; %nearby 2 main (angular distance unit)
    end
    %normal mode ID
    for k=1:1:input_count      
        VB_ID(i,1) = catalog(i,1);
        VB_ID(i,2) = (rad2deg(acos(stars_in_FOV_ID(1,3)*Main_star_vector(1) + stars_in_FOV_ID(1,4)*Main_star_vector(2) + stars_in_FOV_ID(1,5)*Main_star_vector(3)))*3600)*16;  % distance 2 first (angular distance unit)
        VB_ID(i,3) = input_count;
        VB_ID(i,k + 3) = catalog(stars_in_FOV_ID(k,1),1);
    end
    %==================miss first nearby star mode========================%
    
    %SPD
    %miss first nearby star mode
    for k=1:1:(input_count-1)-1
        VB_2(i,1) = catalog(i,1);
        VB_2(i,2) = (rad2deg(acos(stars_in_FOV_ID(2,3)*Main_star_vector(1) + stars_in_FOV_ID(2,4)*Main_star_vector(2) + stars_in_FOV_ID(2,5)*Main_star_vector(3)))*3600)*16; % distance to second (angular distance unit)
        VB_2(i,3) = input_count-1;
        VB_2(i,(2*k) + 2) = (rad2deg(acos(stars_in_FOV_ID((k+2),3)*stars_in_FOV_ID(2,3) + stars_in_FOV_ID((k+2),4)*stars_in_FOV_ID(2,4) + stars_in_FOV_ID((k+2),5)*stars_in_FOV_ID(2,5)))*3600)*16; %angular distance unit
        VB_2(i,(2*k) + 3) = (rad2deg(acos(stars_in_FOV_ID((k+2),3)*Main_star_vector(1) + stars_in_FOV_ID((k+2),4)*Main_star_vector(2) + stars_in_FOV_ID((k+2),5)*Main_star_vector(3)))*3600)*16;    %angular distance unit
    end

    %miss first nearby star mode ID
    for k=1:1:input_count-1
        VB_2_ID(i,1) = catalog(i,1);
        VB_2_ID(i,2) = (rad2deg(acos(stars_in_FOV_ID(2,3)*Main_star_vector(1) + stars_in_FOV_ID(2,4)*Main_star_vector(2) + stars_in_FOV_ID(2,5)*Main_star_vector(3)))*3600)*16; % distance to second (angular distance unit)
        VB_2_ID(i,3) = input_count-1;
        VB_2_ID(i,k + 3) = catalog(stars_in_FOV_ID(k+1,1));
    end
end
VB = bubble_sort(VB,length(VB),2,1);
VB_ID = bubble_sort(VB_ID,length(VB_ID),2,1);
VB_2 = bubble_sort(VB_2,length(VB_2),2,1);
VB_2_ID = bubble_sort(VB_2_ID,length(VB_2_ID),2,1);

%normal mode
SPD_fd = fopen('guide_catalogue_1.txt','w+');
for i=1:1:length(VB)
    for j=1:1:length(VB(1,:))
        if VB(i,j) ~= 0
            fprintf(SPD_fd,'%s ', num2str(round(VB(i,j))));
        else
            fprintf(SPD_fd,'%s ', num2str(round(VB(i,j))));
        end
    end
    fprintf(SPD_fd,'\n');
end
fclose(SPD_fd);

SPD_ID_fd = fopen('guide_catalogue_ID_1.txt','w+');
for i=1:1:length(VB_ID)
    for j=1:1:length(VB_ID(1,:))
        if VB_ID(i,j) ~= 0
            fprintf(SPD_ID_fd,'%s ', num2str(round(VB_ID(i,j))));
        else
            fprintf(SPD_ID_fd,'%s ', num2str(round(VB_ID(i,j))));
        end
    end
    fprintf(SPD_ID_fd,'\n');
end
fclose(SPD_ID_fd);
%miss first nearby star mode
SPD_2_fd = fopen('guide_catalogue_2.txt','w+');
for i=1:1:length(VB_2)
    for j=1:1:length(VB_2(1,:))
        if VB_2(i,j) ~= 0
            fprintf(SPD_2_fd,'%s ', num2str(round(VB_2(i,j))));
        else
            fprintf(SPD_2_fd,'%s ', num2str(round(VB_2(i,j))));
        end
    end
    fprintf(SPD_2_fd,'\n');
end
fclose(SPD_2_fd);

SPD_2_ID_fd = fopen('guide_catalogue_ID_2.txt','w+');
for i=1:1:length(VB_2_ID)
    for j=1:1:length(VB_2_ID(1,:)) 
        if VB_2_ID(i,j) ~= 0
            fprintf(SPD_2_ID_fd,'%s ', num2str(round(VB_2_ID(i,j))));
        else
            fprintf(SPD_2_ID_fd,'%s ', num2str(round(VB_2_ID(i,j))));
        end
    end
    fprintf(SPD_2_ID_fd,'\n');
end
fclose(SPD_2_ID_fd);

