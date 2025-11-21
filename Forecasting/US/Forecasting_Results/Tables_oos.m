

%% 'Industrial Production'
disp('------------------ | U.S monthly data: Industrial Production MSPE | ----------------------------')
R3 = zeros(h,1);

for j = SET
R3 = [R3 squeeze(mean(MSPE(:,1,:,j)))./squeeze(mean(MSPE(:,1,:,benchmark)))];
end

R3 = R3(horizon,2:end);R3 = round(R3, 3);  R3(:,5)=[];
T3 = array2table([horizon' R3], 'VariableNames', [{'h'}, modelNames]);disp(T3)

%% 'PCEPI'
disp('--------------------------- | U.S monthly data: PCEPI MSPE | ----------------------------')
R4 = zeros(h,1);

for j=SET
R4 = [R4 squeeze(mean(MSPE(:,2,:,j)))./squeeze(mean(MSPE(:,2,:,benchmark)))];
end
R4 = R4(horizon,2:end);R4 = round(R4, 3);  R4(:,5)=[];
T4 = array2table([horizon' R4], 'VariableNames', [{'h'}, modelNames]);disp(T4)

%%  'Industrial Production'
disp('--------------------------- | U.S monthly data: Industrial Production QScore90 | ----------------------------')
R7=zeros(h,1);
for j=SET
R7 = [R7 squeeze(mean(QScore90(:,1,:,j)))./squeeze(mean(QScore90(:,1,:,benchmark)))];
end
R7 = R7(horizon,2:end);R7 = round(R7, 3);  R7(:,5)=[];
T7 = array2table([horizon' R7], 'VariableNames', [{'h'}, modelNames]);disp(T7)


%%  'PCEPI'
disp('--------------------------- | U.S monthly data: PCEPI QScore90 | ----------------------------')
R8=zeros(h,1);
for j=SET
R8 = [R8 squeeze(mean(QScore90(:,2,:,j)))./squeeze(mean(QScore90(:,2,:,benchmark)))];
end
R8 = R8(horizon,2:end);R8 = round(R8, 3);  R8(:,5)=[];
T8 = array2table([horizon' R8], 'VariableNames', [{'h'}, modelNames]);disp(T8)

%%  'Industrial Production'
disp('--------------------------- | U.S monthly data: Industrial Production QScore10 | ----------------------------')
R9=zeros(h,1);
for j=SET
R9 = [R9 squeeze(mean(QScore10(:,1,:,j)))./squeeze(mean(QScore10(:,1,:,benchmark)))];
end
R9 = R9(horizon,2:end);R9 = round(R9, 3);  R9(:,5)=[];
T9 = array2table([horizon' R9], 'VariableNames', [{'h'}, modelNames]);disp(T9)


%%  'PCEPI'
disp('--------------------------- | U.S monthly data: PCEPI QScore10 | ----------------------------')
R10=zeros(h,1);
for j=SET
R10 = [R10 squeeze(mean(QScore10(:,2,:,j)))./squeeze(mean(QScore10(:,2,:,benchmark)))];
end
R10 = R10(horizon,2:end);R10 = round(R10, 3);  R10(:,5)=[];
T10 = array2table([horizon' R10], 'VariableNames', [{'h'}, modelNames]);disp(T10)










