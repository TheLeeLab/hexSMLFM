function [model,stdx,mse] = backwardModel(locs,rhoScaling)
% backwardModel

u =locs(:,2);
v = locs(:,3);
x = locs(:,4)-u/rhoScaling;
y = locs(:,5)-v/rhoScaling;
alpha = locs(:,11:12);

b = [x;y];
zeros_t = zeros(size(u,1),1);
ones_t = ones(size(u,1),1);

A=[[ones_t zeros_t alpha(:,1)]; [zeros_t ones_t alpha(:,2)]];

[model,stdx,mse] = lscov(A,b);
