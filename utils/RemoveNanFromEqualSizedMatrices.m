function [A_nanless, B_nanless] = RemoveNanFromEqualSizedMatrices(A,B)
% 
%   UNFINISHED!!! 
% 
% This function removes null entries (nans) from equal sized matrices
%   returning equal sized matrices
% 
% INPUT
%   
%   A            matrix A
%   B            matrix B
% 
% OUTPUT
% 
%   A_nanless    matrix A without nan indices from A and B
%   B_nanless    matrix B without nan indices from A and B
% 
% By Dima Bryzgalov, MOBS team, Paris,
% 07/07/2020
% github.com/bryzgalovdm

%% Arguments handling

% Check number of arguments
if nargin ~= 2
    error('Input should be two matrices');
end

% Check if both matrices are 2D
if length(size(A)) ~= 2 || length(size(B)) ~= 2
    error('Matrices should be two-dimenstional');
end

% Check if matrices A and B are of equal size
if sum(size(A) == size(B)) ~= 2
    error('Matrices should be of equal size');
end

%% Remove nans
idxnan_A = find(sum(A)==0);
idxnan_B = find(sum(B)==0);

end