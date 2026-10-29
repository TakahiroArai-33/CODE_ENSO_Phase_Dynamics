% RECONSTRUCT THE LIFECYCLE OF THE EL NINO SOUTHERN OSCILLATION (ENSO) 
% USING DATA-DRIVEN SPECTRAL ANALYSIS OF KOOPMAN/TRANSFER OPERATORS
%
% Modified 2021/09/14

%% SCRIPT EXECUTION OPTIONS
ifDataSrc           = true; % source data (Indo-Pacific SST)
ifDataTrg           = true; % target data (Nino 3.4 index)
ifNLSA              = true; % perform NLSA (kernel eigenfunctions)
ifKoopman           = true; % compute Koopman eigenfunctions
ifNinoLifecycle     = true;  % ENSO lifecycle from Nino 3.4 index 
ifKoopmanLifecycle  = true;  % ENSO lifecycle from generator eigenfunction
% ifPlotSpectrum      = true;  % plot spectrum of the generator
% ifPlotEnsoLifecycle = true;
% ifPrintFig          = true;  % print figures to file

%% NLSA PARAMETERS
% ERSSTv5 reanalysis
NLSA.dataset           = 'ersstV5';
NLSA.period            = {'197001' '202512'};
NLSA.climatologyPeriod = {'198101' '201012'};
NLSA.sourceVar         = {'IPSST'};
NLSA.targetVar         = {'Nino3.4'};
NLSA.embWindow         = 48; 
NLSA.kernel            = 'cone';
NLSA.ifDen             = false;

experiment = demoKoopman_experiment(NLSA);

%% BATCH PROCESSING
iProc = 1; % index of batch process for this script
nProc = 1; % number of batch processes


%% FIGURE DIRECTORY
figDir = fullfile(pwd, 'figs', experiment);
if ~isdir(figDir)
    mkdir(figDir)
end

%% EXTRACT SOURCE DATA
if ifDataSrc
    for iVar = 1 : numel(NLSA.sourceVar)
        msgStr = sprintf('Reading source data %s...', ...
                          NLSA.sourceVar{iVar});
        disp(msgStr) 
        t = tic;
        demoKoopman_data(NLSA.dataset, NLSA.period, NLSA.sourceVar{iVar}, ...
                          NLSA.climatologyPeriod) 
        toc(t)
    end
end

%% EXTRACT TARGET DATA
if ifDataTrg
    for iVar = 1 : numel(NLSA.targetVar)
        msgStr = sprintf('Reading target data %s...', ...
                           NLSA.targetVar{iVar});
        disp(msgStr) 
        t = tic;
        demoKoopman_data(NLSA.dataset, NLSA.period, NLSA.targetVar{iVar}, ...
                          NLSA.climatologyPeriod) 
        toc(t)
    end
end


%% BUILD NLSA MODEL, DETERMINE BASIC ARRAY SIZES
% In is a data structures containing the NLSA parameters. 
%
% nSE is the number of training samples available after Takens delay
% embedding.
%
% nE is equal to half of the delay embedding window length employed in the 
% training data. 
%
% nSB is the number of samples left out in the start of the time interval (for
% temporal finite differnences employed in the kernel).
%
% nShiftTakens is the temporal shift applied to align Nino indices with the
% center of the Takens embedding window eployed in the NLSA kernel. 

[model, In] = demoKoopman_nlsaModel(NLSA);

nSE = getNTotalEmbSample(model);    
nSB = getNXB(model.embComponent);
nShiftTakens = floor(getEmbeddingWindow(model.embComponent) / 2);

%% SPECIFY GLOBAL PARAMETERS
% The following variables are defined:
% nShiftNino:   Temporal shift to obtain 2D Nino index
% idxZEnso:     ENSO eigenfunction from generator      
% Spec:         Parameters for generator spectral plots
% iC<VarName>:  Index of variable VarName in the target data of the nlsaModel 
% iRec<Name>:   Index of Koopman-reconstructed mode Name in nlsaModel 

nShiftNino = 11;        
iCNino34   = 1;  % Nino 3.4 index

% Plot color for eigenvalue groups
Spec.Color.const   = [0 0 0];
Spec.Color.annual  = [0 0 1];
Spec.Color.semi    = [0.9290 0.6940 0.1250];
Spec.Color.tri     = [1 0 1];
Spec.Color.trend   = [0 1 0];
Spec.Color.trendC  = [0 1 1];
Spec.Color.enso    = [1 0 0];
Spec.Color.ensoC   = [0.4940 0.1840 0.5560]; 
Spec.Color.decadal = [1 1 0] * .5;

Spec.legend = {'seasonal branch' ...
                'trend branch' ...
                'ENSO branch' ...
                'constant' ...
                'annual' ...
                'semiannual' ...
                'triannual' ...
                'trend' ... 
                'trend combination' ...
                'fundamental ENSO' ...
                'ENSO combination' ...
                'decadal' ...
              };


switch experiment

case 'ersstV5_197001-202512_IPSST_emb48_cone'

    datasetStr = 'ERSSTv5';

    % ENSO eigenfunction index
    %* ENSO modeのeigenfunctionの番号
    % 2025-12 まで延ばすと ENSO 周波数の mode は g9/g10 に移る。
    idxZEnso     = 9;


    % Manual spectrum-label indices are disabled in MyCodeOutputData.m.
    % This script only saves Z/Nino34/Lambda; mode labels are checked elsewhere.
    % Spec.Idx.const        = 1;
    % Spec.Idx.annual       = [2 3];
    % Spec.Idx.semi         = [4 5];
    % Spec.Idx.tri          = [13 14];
    % Spec.Idx.trend        = 6;
    % Spec.Idx.trendC       = [9 10];
    % Spec.Idx.enso         = [9 10];
    % Spec.Idx.ensoC        = [11 12 16 17];
    % Spec.Idx.decadal      = 15;
    % Spec.Idx.annualBranch = [13 4 2 3 5 14];
    % Spec.Idx.trendBranch  = [9 6 10];
    % Spec.Idx.ensoBranch   = [16 11 7 8 12 17];

    Spec.xLim  = [-1 .1];
    Spec.yLim  = [-3 3]; 


case 'ersstV4_197001-202002_IPSST_emb48_cone'

    datasetStr = 'ERSSTv4'; 

    % ENSO eigenfunction index
    %* ENSO modeのeigenfunctionの番号
    idxZEnso     = 7;


    % Manual spectrum-label indices are disabled in MyCodeOutputData.m.
    % This script only saves Z/Nino34/Lambda; mode labels are checked elsewhere.
    % Spec.Idx.const        = 1;
    % Spec.Idx.annual       = [2 3];
    % Spec.Idx.semi         = [4 5];
    % Spec.Idx.tri          = [13 14];
    % Spec.Idx.trend        = 6;
    % Spec.Idx.trendC       = [9 10];
    % Spec.Idx.enso         = [7 8];
    % Spec.Idx.ensoC        = [11 12 16 17];
    % Spec.Idx.decadal      = 15;
    % Spec.Idx.annualBranch = [13 4 2 3 5 14];
    % Spec.Idx.trendBranch  = [9 6 10];
    % Spec.Idx.ensoBranch   = [16 11 7 8 12 17];

    Spec.xLim  = [-1 .1];
    Spec.yLim  = [-3 3]; 

case 'hadisst1_197001-202002_IPSST_emb48_cone'

    datasetStr = 'HadISST1';

    % ENSO eigenfunction index
    %* ENSO modeのeigenfunctionの番号
    idxZEnso     = 8;


    % Manual spectrum-label indices are disabled in MyCodeOutputData.m.
    % This script only saves Z/Nino34/Lambda; mode labels are checked elsewhere.
    % Spec.Idx.const        = 1;
    % Spec.Idx.annual       = [2 3];
    % Spec.Idx.semi         = [4 5];
    % Spec.Idx.tri          = [13 14];
    % Spec.Idx.trend        = 6;
    % Spec.Idx.trendC       = [9 10];
    % Spec.Idx.enso         = [7 8];
    % Spec.Idx.ensoC        = [11 12 16 17];
    % Spec.Idx.decadal      = 15;
    % Spec.Idx.annualBranch = [13 4 2 3 5 14];
    % Spec.Idx.trendBranch  = [9 6 10];
    % Spec.Idx.ensoBranch   = [16 11 7 8 12 17];

    Spec.xLim  = [-1 .1];
    Spec.yLim  = [-3 3]; 

case 'hadisst1_197001-202512_IPSST_emb48_cone'

    datasetStr = 'HadISST1';

    % ENSO eigenfunction index
    %* ENSO modeのeigenfunctionの番号
    % 1970-2025 の HadISST1 では ENSO 周波数の mode は g9/g10。
    idxZEnso     = 9;


    % Manual spectrum-label indices are disabled in MyCodeOutputData.m.
    % This script only saves Z/Nino34/Lambda; mode labels are checked elsewhere.
    % Spec.Idx.const        = 1;
    % Spec.Idx.annual       = [2 3];
    % Spec.Idx.semi         = [4 5];
    % Spec.Idx.tri          = [13 14];
    % Spec.Idx.trend        = 8;
    % Spec.Idx.trendC       = [];
    % Spec.Idx.enso         = [9 10];
    % Spec.Idx.ensoC        = [18 19];
    % Spec.Idx.decadal      = 15;
    % Spec.Idx.annualBranch = [13 4 2 3 5 14];
    % Spec.Idx.trendBranch  = [8];
    % Spec.Idx.ensoBranch   = [18 9 10 19];

    Spec.xLim  = [-1 .1];
    Spec.yLim  = [-3 3]; 

case 'cobe2_197001-202002_IPSST_emb48_cone'

    datasetStr = 'COBE-SST2';

    % ENSO eigenfunction index
    %* ENSO modeのeigenfunctionの番号
    idxZEnso     = 9;


    % Manual spectrum-label indices are disabled in MyCodeOutputData.m.
    % This script only saves Z/Nino34/Lambda; mode labels are checked elsewhere.
    % Spec.Idx.const        = 1;
    % Spec.Idx.annual       = [2 3];
    % Spec.Idx.semi         = [4 5];
    % Spec.Idx.tri          = [13 14];
    % Spec.Idx.trend        = 6;
    % Spec.Idx.trendC       = [9 10];
    % Spec.Idx.enso         = [7 8];
    % Spec.Idx.ensoC        = [11 12 16 17];
    % Spec.Idx.decadal      = 15;
    % Spec.Idx.annualBranch = [13 4 2 3 5 14];
    % Spec.Idx.trendBranch  = [9 6 10];
    % Spec.Idx.ensoBranch   = [16 11 7 8 12 17];

    Spec.xLim  = [-1 .1];
    Spec.yLim  = [-3 3]; 

case 'cobe2_197001-202512_IPSST_emb48_cone'

    datasetStr = 'COBE-SST2';

    % ENSO eigenfunction index
    %* ENSO modeのeigenfunctionの番号
    idxZEnso     = 11;


    % Manual spectrum-label indices are disabled in MyCodeOutputData.m.
    % This script only saves Z/Nino34/Lambda; mode labels are checked elsewhere.
    % Spec.Idx.const        = 1;
    % Spec.Idx.annual       = [2 3];
    % Spec.Idx.semi         = [4 5];
    % Spec.Idx.tri          = [13 14];
    % Spec.Idx.trend        = 6;
    % Spec.Idx.trendC       = [9 10];
    % Spec.Idx.enso         = [7 8];
    % Spec.Idx.ensoC        = [11 12 16 17];
    % Spec.Idx.decadal      = 15;
    % Spec.Idx.annualBranch = [13 4 2 3 5 14];
    % Spec.Idx.trendBranch  = [9 6 10];
    % Spec.Idx.ensoBranch   = [16 11 7 8 12 17];

    Spec.xLim  = [-1 .1];
    Spec.yLim  = [-3 3]; 

case 'ersstV5_197001-202002_IPSST_emb48_cone'

    datasetStr = 'ERSSTv5';

    % ENSO eigenfunction index
    %* ENSO modeのeigenfunctionの番号
    idxZEnso     = 7;


    % Manual spectrum-label indices are disabled in MyCodeOutputData.m.
    % This script only saves Z/Nino34/Lambda; mode labels are checked elsewhere.
    % Spec.Idx.const        = 1;
    % Spec.Idx.annual       = [2 3];
    % Spec.Idx.semi         = [4 5];
    % Spec.Idx.tri          = [13 14];
    % Spec.Idx.trend        = 6;
    % Spec.Idx.trendC       = [9 10];
    % Spec.Idx.enso         = [7 8];
    % Spec.Idx.ensoC        = [11 12 16 17];
    % Spec.Idx.decadal      = 15;
    % Spec.Idx.annualBranch = [13 4 2 3 5 14];
    % Spec.Idx.trendBranch  = [9 6 10];
    % Spec.Idx.ensoBranch   = [16 11 7 8 12 17];

    Spec.xLim  = [-1 .1];
    Spec.yLim  = [-3 3]; 

otherwise
    error('Invalid experiment')

end
    
%% PERFORM NLSA 
% Output from each step is saved on disk.
if ifNLSA

    disp('Takens delay embedding for source data...'); t = tic; 
    computeDelayEmbedding(model)
    toc(t)

    % The following step is needed only if we are using velocity-dependent
    % kernels.
    if isa(model.embComponent, 'nlsaEmbeddedComponent_xi')
        disp('Phase space velocity (time tendency of data)...'); t = tic; 
        computeVelocity(model)
        toc(t)
    end

    % The following steps are needed only if we are using variable-bandwidth
    % kernels.
    if isa(model, 'nlsaModel_den')
        fprintf('Pairwise distances for density data, %i/%i...\n', ...
                  iProc, nProc); 
        t = tic;
        computeDenPairwiseDistances(model, iProc, nProc)
        toc(t)

        disp('Distance normalization for kernel density steimation...');
        t = tic;
        computeDenBandwidthNormalization(model);
        toc(t)

        disp('Kernel bandwidth tuning for density estimation...'); t = tic;
        computeDenKernelDoubleSum(model);
        toc(t)

        disp('Kernel density estimation...'); t = tic;
        computeDensity(model);
        toc(t)

        disp('Takens delay embedding for density data...'); t = tic;
        computeDensityDelayEmbedding(model);
        toc(t)
    end

    fprintf('Pairwise distances (%i/%i)...\n', iProc, nProc); t = tic;
    computePairwiseDistances(model, iProc, nProc)
    toc(t)

    disp('Distance symmetrization...'); t = tic;
    symmetrizeDistances(model)
    toc(t)

    disp('Kernel bandwidth tuning...'); t = tic;
    computeKernelDoubleSum(model)
    toc(t)

    disp('Kernel eigenfunctions...'); t = tic;
    computeDiffusionEigenfunctions(model)
    toc(t)
end


%% COMPUTE EIGENFUNCTIONS OF KOOPMAN GENERATOR
if ifKoopman
    disp('Koopman eigenfunctions...'); t = tic;
    computeKoopmanEigenfunctions(model, 'ifLeftEigenfunctions', true)
    toc(t)
end

%% CONSTRUCT NINO-BASED ENSO LIFECYCLE
% Build a structure Nino34 such that:
% 
% Nino34.idx is an array of size [2 nSE], where nSE is the number of samples 
% after delay embedding. Nino34.idx(1, :) contains the values of the 
% Nino 3.4 index at the current time. Nino34(2, :) contains the values of 
% the Nino 3.4 index at nShiftNino timesteps (months) in the past.
% 
% Nino34.time is an array of size [1 nSE] containing the timestamps in
% Matlab serial date number format. 
if ifNinoLifecycle

    disp('Constructing Nino 3.4-based bivariate ENSO index...'); t = tic;

    % Timestamps
    Nino34.time = getTrgTime(model); 
    Nino34.time = Nino34.time(nSB + 1 + nShiftTakens : end);
    Nino34.time = Nino34.time(1 : nSE);

    % Nino 3.4 index
    nino = getData(model.trgComponent(iCNino34));
    Nino34.idx = [nino(nShiftNino + 1 : end) 
                 nino(1 : end - nShiftNino)];
    Nino34.idx = Nino34.idx(:, nSB + nShiftTakens - nShiftNino + 1 : end);
    Nino34.idx = Nino34.idx(:, 1 : nSE);
    Nino34.idx = complex(Nino34.idx(1, :), Nino34.idx(2, :));
    Nino34.angle = angle(Nino34.idx);
end


%! Nino出力
mkdir('myFolder');
save('./myFolder/Nino34_data.mat', 'Nino34');




%% ENSO LIFECYCLE BASED ON KOOPMAN EIGENFUNCTIONS
% Build a structure Z such that:
% 
% Z.idx is a complex array of size [2 nSE], where nSE is the number of 
% samples after delay embedding. Z.idx(1, :) and Z.idx(2, :) contain the 
% real and imaginary parts of the generator-based ENSO index. 
% 
% Z.time is an array of size [1 nSE] containing the timestamps in
% Matlab serial date number format. 

if ifKoopmanLifecycle
    %* Koopman eigenfunctionを持つZを作っている!
    %* Z.idxがidxZEnsoで指定したモードのダイナミクスになる。
    disp('Constructing Koopman-based complex ENSO index...'); t = tic;

    % Retrieve Koopman eigenfunctions
    z = getKoopmanEigenfunctions(model); %* right Koopman eigenfunction
    [z_left, mu] = getLeftKoopmanEigenfunctions(model); %* left Koopman eigenfunction and weights
    T = getEigenperiods(model.koopmanOp); %* 固有周期
    TEnso = abs(T(idxZEnso) / 12); %* periodを年に変換
    a = real(z(:, :))';
    b = imag(z(:, :))';
    Z.idx = complex(a, b);
    Z.time = getTrgTime(model);
    Z.time = Z.time(nSB + 1 + nShiftTakens : end);
    Z.time = Z.time(1 : nSE);
    Z.angle = angle(Z.idx);

    Z_left.idx = complex(real(z_left(:, :))', imag(z_left(:, :))');
    Z_left.time = Z.time;
    Z_left.angle = angle(Z_left.idx);
end


%! eigenfunction出力
mkdir('myFolder');
save('./myFolder/Z.mat', 'Z');
save('./myFolder/Z_left.mat', 'Z_left');
save('./myFolder/mu.mat', 'mu');



% Retrieve Koopman eigenvalues. Set frequency units to cycles/year.
%* 固有値
    Lambda.value = getKoopmanEigenvalues(model);
    Lambda.gamma = Lambda.value * 12 / 2 / pi; 
%! eigenvalue出力
mkdir('myFolder');
save('./myFolder/Lambda.mat', 'Lambda');
