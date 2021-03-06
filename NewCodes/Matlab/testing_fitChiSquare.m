clear all; close all; clc;

xData = 1:10;
yData = [2 2 3 4 7 6 3 2 2 2];

addpath('S:\MC_IDS\Matlab Code\Core Codes');

modelFun = @singletFun;

guess = [5 5 1 1.5];

dxData = [];
dyData = [];
options.Plot = 1;
options.ErrorsUnknown = 1;
options.Display = 'off';

[params,dParams,gof,stddev] = fitChiSquare(xData,yData,modelFun,guess,dxData,dyData,options);
