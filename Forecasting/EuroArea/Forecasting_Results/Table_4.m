%% Replication Table 7
load("oos_EA_p2.mat")
modelNames = {'AVP-VAR','CP-VAR','CP-VAR SV','TVP-VAR-EB','TVP-VAR','VAR SVO-t','FAVAR','FAVAR SV'};
SET = [1 2 4 3 5 9 6 7 8];
benchmark = 5;

%% Table 3: 'GDPC1'
disp('GDP MSPE')
R3=zeros(h,1);
for j=SET
R3 = [R3 squeeze(mean(MSPE(:,1,:,j)))./squeeze(mean(MSPE(:,1,:,benchmark)))];
end
R3 = R3(:,2:end);R3 = round(R3, 3);  R3(:, 5) = [];
T3 = array2table(R3, 'VariableNames', modelNames);disp(T3)

%% Table 4: 'HICP'
disp('HICP MSPE')
R4=zeros(h,1);
for j=SET
R4 = [R4 squeeze(mean(MSPE(:,2,:,j)))./squeeze(mean(MSPE(:,2,:,benchmark)))];
end
R4 = R4(:,2:end);R4 = round(R4, 3);  R4(:, 5) = [];
T4 = array2table(R4, 'VariableNames', modelNames);disp(T4)

%% Table 11: computational cost
SET = [1 2 4 3 9 6 7 8];
disp('Computational cost')
T11 = array2table([sum(elapsedTime(:,SET))./sum(elapsedTime(:,3)) ;max(elapsedTime(:,SET))./60], 'VariableNames', modelNames);disp(T11)






