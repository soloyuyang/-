% PURPOSE: An example of using distance_wm_band
%          Creates a distance based weight matrix
%          based on a band of distances
%---------------------------------------------------
% USAGE: distance_wmd_band
%---------------------------------------------------

clear all
clc
load ct_data    % Load Connecticut county latitude and longitude coordinates
                % 8 counties in all
                % http://en.wikipedia.org/wiki/List_of_counties_in_Connecticut

% The following creates a distance based weight matrix if neighbors are
% between 30 and 60 miles of each other

result = distance_wm_band(yc,xc,30,60);    % Make sure that you enter latitude and longitude
                                           % in that order
                                           

fprintf('Matrix of pairwise distances')
result.dista  % This martrix contains the pairwise distances

fprintf('Lower distance band')
result.dl

fprintf('Upper distance band')
result.du

fprintf('Connectivity Matrix')
result.dw    % This matrix contains the connectivity structure
             % 1 = neighbor, 0 = not a neighbor

fprintf('Sparse row-normalized weight matrix')
result.W     % Sparse row-normalized W matrix 

spy(result.W)



                                