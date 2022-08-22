function HalfLife=GetHalfLife(Isotope)
%
% HalfLife=GetHalfLife(Isotope)
%
% Returns half life for knows Isotopes, for the moment 'C11', 'F18', 'O15',
% 'I123'
%
% HalfLife - half life of substance in minutes
%
% CS, 20140112
%
switch Isotope
    case 'C11'
        HalfLife=20.4;
    case 'F18'
        HalfLife=109.8;
    case 'O15'
        HalfLife=2.05;
    case 'I123'
        HalfLife=792;
    otherwise
        warning('GetHalfLife: Unknown isotope, no correction');
        HalfLife=1e100;
end
