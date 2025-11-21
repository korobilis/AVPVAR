%% Table 2
% Each row is a diferent forecasting horizon
clear; close all; clc;
load("oos_US_p2.mat")
modelNames = {'AVP-VAR','CP-VAR','CP-VAR SV','TVP-VAR-EB','TVP-VAR','VAR SVO-t','FAVAR','FAVAR SV'};
SET = [1 2 4 3 5 9 6 7 8 ];
benchmark = 5; % OLS: benchmark = 5
hh = [1:6 9 12 15 18 24]';


%% Table 3: 'Industrial Production'
disp('Industrial Production MSPE')
R3=zeros(size(hh,1),1);
for i = hh
    for j = SET
        R3 = [R3 squeeze(mean(MSPE(:,1,i,j)))./squeeze(mean(MSPE(:,1,i,benchmark)))];
    end
end
R3 = R3(:,2:end);R3 = round(R3, 3); R3(:,5)=[];  % Round to 3 decimal places
T3 = array2table([hh R3], 'VariableNames', [{'h'},modelNames]);disp(T3)


%% Table 4: 'PCEPI'
disp('PCEPI MSPE')
R4 = zeros(size(hh,1),1);
for i = hh
    for j=SET
        R4 = [R4 squeeze(mean(MSPE(:,2,i,j)))./squeeze(mean(MSPE(:,2,i,benchmark)))];
    end
end
R4 = R4(:,2:end);R4 = round(R4, 3); R4(:,5)=[];  % Round to 3 decimal places
T4 = array2table([hh R4], 'VariableNames', [{'h'},modelNames]);disp(T4)




%% Table 11: computational cost
disp('Computational cost')
SET = [1 2 4 3 9 6 7 8 ];
T11 = array2table([sum(elapsedTime(:,SET))./sum(elapsedTime(:,3)) ;max(elapsedTime(:,SET))./60], 'VariableNames', modelNames);disp(T11)










