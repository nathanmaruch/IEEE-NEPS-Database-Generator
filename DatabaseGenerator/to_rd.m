function [c] = to_rd(m,d)
% -- to_rd(m, d)
% Converts arrays containing the mag (m) and phase (d, in degrees) of
% a complex number to standard complex form
%
% See also: to_pd, to_p, to_r
 if nargin~=2
 error('incorrect number of arguments.');
 end
 % Two arguments: m contains magnitudes, d angles
 c = m .* exp(j*d2r(d));
end
function [m, d] = to_pd(c, N)
% -- to_pd(c, N)
% Returns array(s) containing the mag and phase (in degrees) of
% a complex number.
%
% With two output arguments, the magnitude and phase are returned.
% With one output argument and one input argument the magnitude and
% phase are returned in a hybrid matrix with the magnitude and phase
% adjoined in a single matrix of the form [ MAG : PHASE ].
% With one output argument and two input arguments, the
% magnitude is returned if N == 1 and the phase is returned if N == 2.
%
% See also: to_rd, to_p, to_r
 if (nargout == 2)
 m = abs(c);
 d = r2d(angle(c));
 else
 if (nargin == 2 && N == 1)
 m = abs(c);
 elseif (nargin == 2 && N == 2)
 m = r2d(angle(c))
 else
 m = [abs(c) r2d(angle(c))];
 end
 end
end