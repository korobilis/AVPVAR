function [R5,R6,R7,R8,R9,R10] = Tables_oos_function(loadFile)
modelNames = {'AVP-VAR','CP-VAR','CP-VAR SV','TVP-VAR-EB','OLS-iterative','TVP-VAR','VAR SVO-t','FAVAR','FAVAR SV'};
SET = [1 2 4 3 5 9 6 7 8];
benchmark = 5;
load(loadFile);



%% Table 5: 'GDPC1'
disp('GDPC1 MAE')
R5=zeros(h,1);
for j=SET
R5 = [R5 squeeze(mean(MAE(:,1,:,j)))./squeeze(mean(MAE(:,1,:,benchmark)))];
end
R5 = R5(:,2:end);R5 = round(R5, 3);  % Round to 3 decimal places
T5 = array2table(R5, 'VariableNames', modelNames);disp(T5)
%RR=R5; plot_tables_MAE

%% Table 6: 'PCECTPI'
disp('PCECTPI MAE')
R6=zeros(h,1);
for j=SET
R6 = [R6 squeeze(mean(MAE(:,2,:,j)))./squeeze(mean(MAE(:,2,:,benchmark)))];
end
R6 = R6(:,2:end);R6 = round(R6, 3);  % Round to 3 decimal places
T6 = array2table(R6, 'VariableNames', modelNames);disp(T6)
%RR=R6; plot_tables_MAE

%% Table 7: 'GDPC1'
disp('GDPC1 QScore90')
R7=zeros(h,1);
for j=SET
R7 = [R7 squeeze(mean(QScore90(:,1,:,j)))./squeeze(mean(QScore90(:,1,:,benchmark)))];
end
R7 = R7(:,2:end);R7 = round(R7, 3);  % Round to 3 decimal places
T7 = array2table(R7, 'VariableNames', modelNames);disp(T7)
%RR=R7; plot_tables_QS90

%% Table 8: 'PCECTPI'
disp('PCECTPI QScore90')
R8=zeros(h,1);
for j=SET
R8 = [R8 squeeze(mean(QScore90(:,2,:,j)))./squeeze(mean(QScore90(:,2,:,benchmark)))];
end
R8 = R8(:,2:end);R8 = round(R8, 3);  % Round to 3 decimal places
T8 = array2table(R8, 'VariableNames', modelNames);disp(T8)
% RR=R8; plot_tables_QS90

%% Table 9: 'GDPC1'
disp('GDPC1 QScore10')
R9=zeros(h,1);
for j=SET
R9 = [R9 squeeze(mean(QScore10(:,1,:,j)))./squeeze(mean(QScore10(:,1,:,benchmark)))];
end
R9 = R9(:,2:end);R9 = round(R9, 3);  % Round to 3 decimal places
T9 = array2table(R9, 'VariableNames', modelNames);disp(T9)
% RR=R9; plot_tables_QS10

%% Table 10: 'PCECTPI'
disp('PCECTPI QScore10')
R10=zeros(h,1);
for j=SET
R10 = [R10 squeeze(mean(QScore10(:,2,:,j)))./squeeze(mean(QScore10(:,2,:,benchmark)))];
end
R10 = R10(:,2:end);R10 = round(R10, 3);  % Round to 3 decimal places
T10 = array2table(R10, 'VariableNames', modelNames);disp(T10)
% RR=R10; plot_tables_QS10

%% Table 11: computational cost
disp('Computational cost')
T11 = array2table([sum(elapsedTime(:,SET))./sum(elapsedTime(:,3)) ;max(elapsedTime(:,SET))./60], 'VariableNames', modelNames);disp(T11)

R5(:, 5) = [];R6(:, 5) = [];R7(:, 5) = [];R8(:, 5) = [];R9(:, 5) = [];R10(:, 5) = [];


end


