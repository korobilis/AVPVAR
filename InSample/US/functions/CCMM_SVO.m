function [ydraws]   = CCMM_SVO(Y,sample,ydates,p,yrealized,nfore)
tic;
%% Settings for CCMM VAR-SVO code
N = size(Y,2);
ndxYIELDS = [];
ELBbound = [];
check_stationarity = 0;
minnesotaPriorMean = zeros(N,1);
doTightPrior = false;
doRobustPrior = false;
np = 12;
fcstNdraws = 5000;   % 10000
MCMCdraws  = 1000; % 1000
fcstNhorizons = nfore;
rndStream = RandStream.getGlobalStream();

%% SVO prior
SVOalpha = 1 / ( 4 * np) * 10 * np;
SVObeta  = 10 * np - SVOalpha;
SVOmaxscale         = 20;

[PAI_all, PHI_all, invA_all, sqrtht_all, ...
        SVOprob_all, SVOscale_all, ...
        ydraws, yhat, yhatRB, fcstSVdraws, fcstSVoutliers, logscoredraws, logscoredrawsRB] ...
        = mcmcVARSVO(sample, MCMCdraws, p, np, Y, ydates, ...
        minnesotaPriorMean, doTightPrior, doRobustPrior, ...
        SVOalpha, SVObeta, SVOmaxscale, ...
        ndxYIELDS, ELBbound, ...
        check_stationarity, ...
        yrealized,...
        fcstNdraws, fcstNhorizons, rndStream); 
toc