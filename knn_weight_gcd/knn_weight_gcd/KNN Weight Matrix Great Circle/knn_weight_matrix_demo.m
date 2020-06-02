%% Program to demonstrate the knn_weight_matrix code
% Uses knn_weight_matrix.m and haversine.m code

%% Load demonstration data
% Latitude/Longitude coordiantes of Connecticut counties
clear 
clc

load ct_data

% lat = latitude coordinates
% long = longitude coordinates

%% Build 2 nearest neighbor weight matrix

% W is the sparse representation
W = knn_weight_matrix_gcd(xc,yc,2)

% WW is the full (non-sparse) version
WW = full(W)
