function[correct,num_matched,final_output,num_mismatch,final_mainstar] = TRK_matching(num_input_stars,mainstar,nearby_star_vector,searching_catalogue,searching_catalogue_ID,error_range,limit_matching_number)

        num_matched=0;
        correct=1;

     if limit_matching_number > num_input_stars - 1
         limit_matching_number = num_input_stars - 1;
     end
     
    
    % distance by first nearby star to mainstar
    dis_first2main = acos(nearby_star_vector(1,1)*mainstar(1) + nearby_star_vector(1,2)*mainstar(2) + nearby_star_vector(1,3)*mainstar(3));
    
    %SPD in angle
    for k=1:1:(num_input_stars-2)
        VB(1,(2*k) - 1) = (rad2deg(acos(nearby_star_vector(k+1,1)*nearby_star_vector(1,1) + nearby_star_vector(k+1,2)*nearby_star_vector(1,2) + nearby_star_vector(k+1,3)*nearby_star_vector(1,3)))*3600)*16;
        VB(1,(2*k)) = (rad2deg(acos(nearby_star_vector(k+1,1)*mainstar(1) + nearby_star_vector(k+1,2)*mainstar(2) + nearby_star_vector(k+1,3)*mainstar(3)))*3600)*16;
    end
    
        %======================matching=======================%
        
        %=======vote==========%
    vote_list = zeros(length(searching_catalogue(:,1)),num_input_stars); %vote list
    search_count = 0;
    
    for i = 1:1:length(searching_catalogue(:,1))
        search_count = search_count + 1;
        catalogue_search_count = 0;

        for j =1:1:(num_input_stars-2)
            catalogue_search_count = catalogue_search_count + 1;
            if catalogue_search_count > searching_catalogue(i,3)-1
                catalogue_search_count = catalogue_search_count - 1;
                break;
            end
            if searching_catalogue(i,(2*j)+2) == 0 || searching_catalogue(i,(2*j)+3) == 0
                continue;
            end
            ori_flag = 0;
            jump_flag = 0;
                % if route is not match, search next data
            while VB(1,(2*j)-1) < searching_catalogue(i,(2*catalogue_search_count)+2) - error_range || VB(1,(2*j)-1) > searching_catalogue(i,(2*catalogue_search_count)+2) + error_range
                if  ori_flag == 0
                    ori = catalogue_search_count;
                    ori_flag = 1;
                end
                catalogue_search_count = catalogue_search_count + 1;
                if catalogue_search_count > searching_catalogue(i,3)-1
                    jump_flag = 1;
                    catalogue_search_count = ori - 1;
                    break;
                elseif catalogue_search_count > ori + 5
                     jump_flag = 1;
                    catalogue_search_count = ori - 1;
                    break;
                end
            end
            if jump_flag == 1
                continue;
            else
                if VB(1,(2*j)-1) >= searching_catalogue(i,(2*catalogue_search_count)+2) - error_range && VB(1,(2*j)-1) <= searching_catalogue(i,(2*catalogue_search_count)+2) + error_range
                    if VB(1,(2*j)) >= searching_catalogue(i,(2*catalogue_search_count)+3) - error_range && VB(1,(2*j)) <= searching_catalogue(i,(2*catalogue_search_count)+3) + error_range
                        vote_list(search_count,1) = vote_list(search_count,1) + 1;
                        vote_list(search_count,2) = searching_catalogue_ID(search_count ,4);
                        vote_list(search_count,j+2) = searching_catalogue_ID(search_count ,catalogue_search_count +4);
                    end
                end
            end
        end
    end
    %find most votes candidate
    [max_vote_list,index_vote_list] = max(vote_list(:,1));
    
     %deal with same votes candidate
   candidate_number = 0;
    for i=1:1:length(vote_list(:,1))
        if vote_list(i,1) == vote_list(index_vote_list,1)
            candidate_number = candidate_number + 1;
            candidate(candidate_number,:) = [i vote_list(i,:)]; %[vote       ID]
        end
    end
    same_vote_constrain = 10;
    if length(candidate(1,:)) < same_vote_constrain
        same_vote_constrain = length(candidate(1,:));
    end
    if candidate_number > 1
        vote_list_2 = zeros(candidate_number,2);    %[index_vote_list      first 0 position]
        for i=1:1:candidate_number
            vote_list_2(i,1) = candidate(i,1);
            for j=3:1:same_vote_constrain
                if candidate(i,j) == 0
                    vote_list_2(i,2) = vote_list_2(i,2) + 1;    %zero count
                end
            end
        end
        vote_list_2 = bubble_sort(vote_list_2,candidate_number,2,1);
        index_vote_list = vote_list_2(1,1);
    end
    
    % check match results
    final_mainstar = 0;
    num_mismatch = 0;
    if searching_catalogue_ID(index_vote_list,1) ~= mainstar(4)
        final_output = zeros(1,num_input_stars);
        final_mainstar = searching_catalogue_ID(index_vote_list,1);
        correct = 0;
    else
        final_mainstar = searching_catalogue_ID(index_vote_list,1);
        final_output = vote_list(index_vote_list,:);
        num_matched = num_matched + 1;  %mainstar match
        for j=1:1:limit_matching_number
            if final_output(1,j+1) == nearby_star_vector(j,4)
                num_matched = num_matched + 1; %nearby star match
            elseif final_output(1,j+1) ~= nearby_star_vector(j,4) && final_output(1,j+1) ~= 0
                num_mismatch = num_mismatch + 1;    %mismatch
                correct = 0;
            end
        end
    end

end