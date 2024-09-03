function [nearby_star_vector,input_count,ori_reliable_nearby_star_number] = use_weight_arrange_star(use_by_weight_arrange_vector,num_nearby_star,input_stars,input_stars_unit_vector)


%=================caculate smallest distance to nearby star===============%
for i=1:1:num_nearby_star
    clear distance;
    for j=1:1:num_nearby_star
        distance(j,1) = acos(use_by_weight_arrange_vector(i,2)*use_by_weight_arrange_vector(j,2) + use_by_weight_arrange_vector(i,3)*use_by_weight_arrange_vector(j,3) + use_by_weight_arrange_vector(i,4)*use_by_weight_arrange_vector(j,4));
        % exclude itself
        if i == j
            distance(j,1) = 10000;      
        end
    end
    distance = bubble_sort(distance,num_nearby_star,1,1);
    use_by_weight_arrange_vector(i,7) = distance(1,1);  %Smallest distance to nearby star(rad) = SD
end

%========================cacilate nearby stars weight==================%
for i=1:1:num_nearby_star
    stars_in_FOV_CSS(i,1) = i;   %which one in use_by_weight_arrange_vector
    [stars_in_FOV_CSS(i,2),stars_in_FOV_CSS(i,3)] = to_weight(use_by_weight_arrange_vector(i,7),use_by_weight_arrange_vector(i,5));
    % PS : Sequence of (stars_in_FOV_CSS) : 1.SD  2.distance to mainstar
    stars_in_FOV_CSS(i,4) = use_by_weight_arrange_vector(i,6);  %star bright
end

%============remove unreliable nearby star and arrange nearby reliable star==========%
fail_count = 0;    %unreliable nearby star
OK_count = 0;   %reliable nearby star
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
    ori_reliable_nearby_star_number = OK_count;
    OK_CSS_star = bubble_sort(OK_CSS_star,OK_count,3,0);
    first_flag = 0;
    input_count = 0;        % final input star number
    if OK_count > 1
        for k=2:1:OK_count
            if first_flag == 0
                input_count = input_count + 1;
                nearby_star_vector(input_count,1) = use_by_weight_arrange_vector(OK_CSS_star(1,1),2);  %x
                nearby_star_vector(input_count,2) = use_by_weight_arrange_vector(OK_CSS_star(1,1),3);  %y
                nearby_star_vector(input_count,3) = use_by_weight_arrange_vector(OK_CSS_star(1,1),4);  %z
                nearby_star_vector(input_count,4) = input_stars(input_stars_unit_vector(use_by_weight_arrange_vector(OK_CSS_star(1,1),1),4),1);  %ID
                first_flag = 1;
            end
            % if same weight, arranged by bright
            if OK_CSS_star(k,3) + 0.008 >= OK_CSS_star(k-1,3)
                if OK_CSS_star(k,4) < OK_CSS_star(k-1,4)
                    input_count = input_count + 1;
                    nearby_star_vector(input_count,1) = nearby_star_vector(input_count - 1,1);
                    nearby_star_vector(input_count,2) = nearby_star_vector(input_count - 1,2);
                    nearby_star_vector(input_count,3) = nearby_star_vector(input_count - 1,3);
                    nearby_star_vector(input_count,4) = nearby_star_vector(input_count - 1,4);
                    
                    nearby_star_vector(input_count - 1,1) = use_by_weight_arrange_vector(OK_CSS_star(k,1),2);
                    nearby_star_vector(input_count - 1,2) = use_by_weight_arrange_vector(OK_CSS_star(k,1),3);
                    nearby_star_vector(input_count - 1,3) = use_by_weight_arrange_vector(OK_CSS_star(k,1),4);
                    nearby_star_vector(input_count - 1,4) = input_stars(input_stars_unit_vector(use_by_weight_arrange_vector(OK_CSS_star(k,1),1),4),1);  %ID
                else
                    input_count = input_count + 1;
                    nearby_star_vector(input_count,1) = use_by_weight_arrange_vector(OK_CSS_star(k,1),2);  %x
                    nearby_star_vector(input_count,2) = use_by_weight_arrange_vector(OK_CSS_star(k,1),3);  %y
                    nearby_star_vector(input_count,3) = use_by_weight_arrange_vector(OK_CSS_star(k,1),4);  %z
                    nearby_star_vector(input_count,4) = input_stars(input_stars_unit_vector(use_by_weight_arrange_vector(OK_CSS_star(k,1),1),4),1);  %ID
                end
            else
                input_count = input_count + 1;
                nearby_star_vector(input_count,1) = use_by_weight_arrange_vector(OK_CSS_star(k,1),2);  %x
                nearby_star_vector(input_count,2) = use_by_weight_arrange_vector(OK_CSS_star(k,1),3);  %y
                nearby_star_vector(input_count,3) = use_by_weight_arrange_vector(OK_CSS_star(k,1),4);  %z
                nearby_star_vector(input_count,4) = input_stars(input_stars_unit_vector(use_by_weight_arrange_vector(OK_CSS_star(k,1),1),4),1);  %ID
            end
        end
    else
            input_count = input_count + 1;
            nearby_star_vector(input_count,1) = use_by_weight_arrange_vector(OK_CSS_star(1,1),2);  %x
            nearby_star_vector(input_count,2) = use_by_weight_arrange_vector(OK_CSS_star(1,1),3);  %y
            nearby_star_vector(input_count,3) = use_by_weight_arrange_vector(OK_CSS_star(1,1),4);  %z
            nearby_star_vector(input_count,4) = input_stars(input_stars_unit_vector(use_by_weight_arrange_vector(OK_CSS_star(1,1),1),4),1);  %ID

    end
end