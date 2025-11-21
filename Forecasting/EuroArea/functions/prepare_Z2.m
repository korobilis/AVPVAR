function [z,M] = prepare_Z2(Z,p)

%% Inputs
% p is the number of lags in the VAR equation.
%% Outputs
% Z: lagged instruments 


M = size(Z,2); 
z = Z(p:end-1,:); 






% if pZ==0
%     Zfore = repmat(Z(end,:), nfore, 1);
% else
%     Zfore = [Z(end-pZ:end,:); repmat(Z(end,:),nfore-pZ,1)];
% end



end
