function [x,y,z] = to_unit_vector(ra,dec)

    x = cosd(ra)*cosd(dec);
    y= sind(ra)*cosd(dec);
    z = sind(dec);
    
end