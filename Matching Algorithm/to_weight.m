function [W1,W2] = to_weight(sd,dis2main)
    min_star_distance = 0.0017; %rad
    FOV_radius = 8.5;
    best_star_distance_for_attitude = 3;

    if sd > min_star_distance
        W1 = 1;     %not double star
    else
        W1 = 0;     %double star
    end
    
    C2 = (FOV_radius - best_star_distance_for_attitude)^2 + 1;
    
    if sd > min_star_distance
        z = -((best_star_distance_for_attitude - sd)^2) + C2;
    else
        z = 0;
    end
    
    y = (FOV_radius + 1 - dis2main);
    
    
    W2 = y;
    %W2 = y*z;
end