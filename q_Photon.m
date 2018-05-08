% Copyright (C) 2014 Quan Wang <wangq10@rpi.edu>, 
% Signal Analysis and Machine Perception Laboratory, 
% Department of Electrical, Computer, and Systems Engineering, 
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA

function Iq=q_Photon(theta)
%% The luminous intensity distribution of our Vivia 7DR3-RGB fixture (from Renaissance Lighting)
%    theta: angle to normal direction
%    Iq: luminous intensity

%% Photon's q function
a = 0:10:90;
I = [106.5714, 103.4286, 100.5714, 89.7143, 65.4286, 27.7143, 13.5714, ...
    12.2857, 9.5714, 6.1429] / 106.5714;

Iq = interp1(a,I,theta/pi*180,'pchip');
