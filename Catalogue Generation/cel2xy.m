function [x3,y3,z3] = cel2xy(ra,dec,ra0,dec0,phi_r)
%fixed angle rorations from ECI to STR
%intial atittude of STR is aligned to ECI
%Jack 2018/11/18
deg2rad = pi/180;


dec = dec * deg2rad;
ra  = ra * deg2rad;
phi_z = (-1)*ra0*deg2rad;
phi_y = -(90-dec0)*deg2rad;
% phi_r = phi_r* deg2rad; % NCKU' origianl code
phi_r = -phi_r* deg2rad; % phi_r should be inverse (2018/11/18 by Jack)  

x=sin(pi/2-dec).*cos(ra);
y=sin(pi/2-dec).*sin(ra);
z=cos(pi/2-dec);


x1 = x*cos(phi_z)-y*sin(phi_z);
y1 = x*sin(phi_z)+y*cos(phi_z);
z1 = z;

x2 = x1*cos(phi_y)+z1*sin(phi_y);
y2 = y1;
z2 = -x1*sin(phi_y)+z1*cos(phi_y);

x3 = x2*cos(phi_r)-y2*sin(phi_r);
y3 = x2*sin(phi_r)+y2*cos(phi_r);
z3 = z2;

