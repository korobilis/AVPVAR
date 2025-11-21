
%% Appendix_Forecasting_US monthly data: All Tables
clear;clc;
modelNames = {'AVP-VAR','CP-VAR','CP-VAR SV','TVP-VAR-EB','TVP-VAR','VAR SVO-t','FAVAR','FAVAR SV'};
SET = [1 2 4 3 5 9 6 7 8 ];
horizon =[1:6 9 12 15 18 24];
benchmark = 5; % OLS: benchmark = 5

disp('===============================| VAR(1) U.S. Monthly Data: Results relative to OLS benchmark | ===============================')
load("oos_US_p1.mat");
Tables_oos

disp('===============================| VAR(2) U.S. Monthly Data: Results relative to OLS benchmark | ===============================')
load("oos_US_p2.mat");
Tables_oos

disp('===============================| VAR(3) U.S. Monthly Data: Results relative to OLS benchmark | ===============================')
load("oos_US_p3.mat");
Tables_oos

disp('===============================| VAR(4) U.S. Monthly Data: Results relative to OLS benchmark | ===============================')
load("oos_US_p4.mat");
Tables_oos
