% PURPOSE: An example of using invdistance_wm
%          Creates an inverse distance based weight matrix
%          with at least one neighbor
%          from latitude and longitude coordinates
%---------------------------------------------------
% USAGE: invdistance_wmd
%---------------------------------------------------

clear all
load ct_data    % Load Connecticut county latitude and longitude coordinates
                % 8 counties in all
                % http://en.wikipedia.org/wiki/List_of_counties_in_Connecticut

result = invdistance_wm(yc,xc);    % Make sure that you enter latitude and longitude
                                   % in that order

fprintf('Matrix of pairwise distances')
result.dista  % This martrix contains the pairwise distances

fprintf('Inverse Distance Matrix')
result.dw    % This matrix contains the inverse distances
             
spy(result.dw)
