function [Vdphi] = get_Vdphi(c1, c2)
% Assume dvphi is passed through amplifier of 100 
% -> 10mv/degree is transformed to 1v/degree

Vdphi = abs(2*pi*1e-2*(sqrt(sum(c1.^2))-sqrt(sum(c2.^2)))*13.56e6/299792458*180/pi);

end

