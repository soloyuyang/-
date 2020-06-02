function W = knn_weight_matrix_gcd(xc,yc,k)

%% PURPOSE: The function builds a k nearest neighbor N x N row-stochastic spatial weight matrix
% The function returns a sparse matrix W representing the spatial weight
% matrix

%% USAGE: W = knn_weight_matrix(xc,yc,k)
% where the following holds:
%           xc = an (n x 1) vector of longitude coordinates (in degrees)
%           yc = an (n x 1) vector of latitude coordinates (in degrees)
%           k = the number of nearest neighbors

%% NOTES: The code uses the Gret Circle distance formula to create an
% N x N matrix of pairwise distances

%% AUTHOR: Donald J. Lacombe
% Donald J. Lacombe
% Associate Professor
% Department of Personal Financial Planning
% Texas Tech University
% donald.lacombe@ttu.edu
% I greatly appreciate the help of Mark Burkey of North Carolina A&T
% University.
% Special thanks to Stuart McIntyre of the University of Strathclyde for
% help in improving the function.

%% Rearrange coordinates and pre-allocate matrices
n = length(xc);
tempw = zeros(n,n);
nnlist = zeros(n,k);


%% Create pairwise distance matrix using Great Circle distance calculation
temp2 = GeogDistance(xc,yc);

%% Find indices for the nearest neighbors
for i=1:n;
    temp1 = temp2(i,:);
    [~,xind] = sort(temp1);
    nnlist(i,1:k) = xind(1,2:k+1);
end

%% Clear the distance matrix
clear temp2

%% Assign a value of 1 to nearest neighbors
for i=1:n
    tempw(i,nnlist(i,:)) = 1;
end

%% Create row-stochastic and sparse spatial weight matrix
W = sparse(tempw./k);

end


%% Great Circle Distance formula

function dista = GeogDistance(lat,long)
[n1,p1] = size(lat);
[n2,p2] = size(long);

if (n1~=n2)
    error('GeogDistance: numbers of locations must be identical for latitude, longitude');
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

