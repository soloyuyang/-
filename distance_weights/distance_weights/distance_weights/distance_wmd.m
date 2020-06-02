% PURPOSE: An example of using distance_wm
%          Creates a distance based weight matrix
%          with at least one neighbor
%          from latitude and longitude coordinates
%---------------------------------------------------
% USAGE: distance_wmd
%---------------------------------------------------

clear all
load ct_data    % Load Connecticut county latitude and longitude coordinates
                % 8 counties in all
                % http://en.wikipedia.org/wiki/List_of_counties_in_Connecticut

result = distance_wm(yc,xc);    % Make sure that you enter latitude and longitude
                                % in that order

fprintf('Matrix of pairwise distances')
result.dista  % This martrix contains the pairwise distances

fprintf('Connectivity Matrix')
result.dw    % This matrix contains the connectivity structure
             % 1 = neighbor, 0 = not a neighbor

fprintf('Sparse row-normalized weight matrix')
result.W     % Sparse row-normalized W matrix 

spy(result.W)



                                