function results = distance_wm_band(lat,long,dl,du)
% PURPOSE: Creates a distance based spatial weight matrix with bands
%          using latutide and longitude coordinates and 
%          the Great Circle distance formula
% -------------------------------------------------------------------------
% USAGE: result = distance_wm(lat,long,dl,du)
% where: lat  = latitude coordinates (decimal)
%        long = longitude coordinates (decimal)
%        dl = lower distance band (miles or km)
%        du = upper distance band (miles or km)
% -------------------------------------------------------------------------
% RETURNS: a structure
%          results.dista = (n x n) matrix of pairwide distances
%          results.W     = (n x n) sparse row-standardized weight matrix
%          results.dw    = (n x n) matrix of neighbors (0,1)
%          results.dl    = lower band value
%          results.du    = upper band value
% -------------------------------------------------------------------------        
% NOTES: The default measurement for distances is in miles. If you wish to
%        use kilometers, comment out the Earth's radius in miles, line 143
%        and uncomment line 142
%       
%        The distance band in inclusive of the cutoff points
% -------------------------------------------------------------------------
% REFERENCES: This function uses Thomas Pingle's code for calculating the
%             Great Circle distance which is available here:
%  http://www.geog.ucsb.edu/~pingel/matlabCode/files/greatCircleDistance.m
%
% This code is designed to be used with Jim LeSage's Spatial Econometrics Toolbox
% http://www.spatial-econometrics.com/
%
% Special thanks to Tony E. Smith of the University of Pennsylvania for
% providing technical assistance
% http://www.seas.upenn.edu/~tesmith/
%
% Special thanks to Justin Ross of Indiana University for helpful
% discussions
% http://jross08.googlepages.com/
% -------------------------------------------------------------------------
%
% written by:
% Donald J. Lacombe
% Associate Professor
% Department of Personal Financial Planning
% Texas Tech University
% donald.lacombe@ttu.edu

if nargin ~= 4
error('Wrong # of input arguments to distance_wm_band');
end;

dista = GeogDistance(lat,long);      % Create pairwise distance matrix

dw = dista;
[n k] = size(dw);
for i = 1:n
    for j = 1:n
        if (dw(i,j)==0)
            dw(i,j)=0;
        elseif (dw(i,j)>=dl && dw(i,j)<=du); % Create matrix based on band
            dw(i,j)=1;                       % of distances   
        else
            dw(i,j)=0;
        end
    end
end

W = sparse(normw(dw));

results.dista = dista; % The original distances in miles or km
results.W     = W;     % The sparse distance weight matrix
results.dw    = dw;    % Returned for testing code
results.dl    = dl;    % Lower band value
results.du    = du;    % Upper band value
end

% -------------------------------------------------------------------------
% Start of Thomas Pingle's code for calcualting pairwise Great Circle
% distances

function dista = GeogDistance(lat,long)
  [n1,p1] = size(lat);
  [n2,p2] = size(long);

  if (n1~=n2)
    error('  GeogDistance: numbers of locations must be identical for latitude, longitude');
  else
    n = n1;
  end;

  if (p1>1)                                     % Convert latitude minutes, seconds to hundredths
    s = sign(lat(:,1));                         %   of degree
    c = abs(lat(:,1)) + abs(lat(:,2)./60);
    if (p1 > 2)
      c = c + abs(lat(:,3)./3600);
    end;
    lat = s.*c;
  end;

  if (p2>1)                                     % Convert longitude minutes, seconds to hundredths
    s = sign(long(:,1));                        %   of degree
    c = abs(long(:,1)) + abs(long(:,2)./60);
    if (p2 > 2)
      c = c + abs(long(:,3)./3600);
    end;
    long = s.*c;
  end;

  dista = zeros(n,n);                  % Allocate output distance matrix

  for i = 1:(n-1)                     % Cycle thru all possible pairs of localities
    lat1 = lat(i);
    long1 = long(i);

    for j = (i+1):n
      lat2 = lat(j);
      long2 = long(j);
      
      degToRad = 2*pi/360;
      a = lat1*degToRad;
      b = long1*degToRad;
      c = lat2*degToRad;
      d = long2*degToRad;
      
      if (abs(a-c)<eps && abs(b-d)<eps)
        d = 0;
      else
        t = sin(a)*sin(c)+cos(a)*cos(c)*cos(b-d);
        %r = 6366.2;                               % Earth's radius in km
         r = 3956.545;                              % Earth's radius in miles
        if (t>1)
          d = r*acos(1);
        else
          d = r*acos(t);
        end;
      end;

      dista(i,j) = d;                     
      dista(j,i) = d;
    end;
  end;

end
