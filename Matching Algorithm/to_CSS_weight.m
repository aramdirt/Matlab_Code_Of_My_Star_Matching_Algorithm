function [W1,W2] = to_CSS_weight(sd,dis2main,bright)
    min_star_distance = 0.0017; %rad
    FOV_length = 8.5;
    best_star_distance_for_attitude = 3;

    if sd > min_star_distance
        W1 = 1;     %not double star
    else
        W1 = 0;     %double star
    end
    
    if sd > min_star_distance
        z = -((best_star_distance_for_attitude - sd)^2) + 21.25;
    else
        z = 0;
    end
    
    y = (FOV_length-dis2main);
    
    
    W2 = y;
    %W2 = y*z;
end