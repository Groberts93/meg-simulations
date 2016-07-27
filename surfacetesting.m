%George Roberts - 27/07/2016
%Surface function testing
clear;

addpath('vfunc');
addpath('sfunc');

load('data/spherepts.mat');

xpts = EEGPts1(:,1);
ypts = EEGPts1(:,2);
zpts = EEGPts1(:,3);
T = delaunay(xpts, ypts, zpts);

trisurf(T, xpts, ypts, zpts);
hold on;


for n = 1:length(EEGPts3)
plot3(EEGPts3(n,1), EEGPts3(n,2), EEGPts3(n,3), 'bx', 'LineWidth', 4);
hold on;
end