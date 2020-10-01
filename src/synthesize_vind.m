function [vind] = synthesize_vind(center,norm)
f=13.56e6;    Z01=50;   Z02=1e6;

reader_coil = COIL( 15e-2, 1, 1);
reader_abcd = RESISTOR('S', 50).ABCD(f);

tag_coil = COIL( 5e-2, 1, 10);
tag_abcd = RESISTOR('S', 1).ABCD(f) * ...
    CAPACITOR('P', 1e-9).ABCD(f) * ...
    RESISTOR('P', 1e3).ABCD(f);


tag_coil.move(center, norm);
M = reader_coil.mutualInductance(tag_coil, f);
ms_abcd = INDUCTOR('S',-M).ABCD(f);
mp_abcd = INDUCTOR('P', M).ABCD(f);
abcd = reader_abcd * ms_abcd * mp_abcd * ms_abcd * tag_abcd;
A = abcd(1,1); B = abcd(1,2); C = abcd(2,1); D = abcd(2,2);
S21 = (2*sqrt(real(Z01)*real(Z02)))./(A*Z02+B+C*Z01*Z02+D*Z01);
vind = abs(S21) * 5e4 - .3;

end

