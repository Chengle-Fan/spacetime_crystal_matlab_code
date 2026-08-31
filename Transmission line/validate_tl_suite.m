function summary = validate_tl_suite()
%VALIDATE_TL_SUITE Regression tests for the transmission-line module.
% Run from this directory with MATLAB R2020a+ and Base MATLAB only.

fprintf('Transmission-line validation started with MATLAB %s.\n',version);
moduleDirectory = fileparts(mfilename('fullpath'));
originalDirectory = pwd;
cleanupDirectory = onCleanup(@() cd(originalDirectory));
cd(moduleDirectory);
assert(startsWith(which('tl_build_model'),moduleDirectory), ...
    'tl_build_model does not resolve to this module directory.');

testNames = { ...
    'paper anchors and assumptions', ...
    'static analytic dispersion and TMM', ...
    'PWE versus time TMM', ...
    'passive loss sign', ...
    'finite-chain FDTD energy and Q continuity', ...
    'CROW finite-chain smoke test', ...
    'synthetic signed x-t FFT', ...
    'synthetic fixed-probe FFT and Floquet folding', ...
    'finite-chain Gaussian k scan', ...
    'Bloch FFT versus TMM', ...
    'fit and component interfaces', ...
    'Code Analyzer'};
passed = false(size(testNames));

%% Paper anchors and explicit assumptions
sspp = tl_build_model(struct('topology','sspp','allowAssumptions',true));
expectedLeffective = 4/((2*pi*371e6)^2*19.8e-12);
assert_relative(sspp.derived.effectiveBoundaryInductance, ...
    expectedLeffective,1e-12,'SSPP paper anchor');
crow = tl_build_model(struct('topology','crow','allowAssumptions',true));
assert_relative(crow.derived.fcolHz,314.7e6,2e-3,'CROW fcol anchor');
assumptionRejected = false;
try
    tl_build_model(struct('topology','sspp'));
catch
    assumptionRejected = true;
end
assert(assumptionRejected,'Missing paper parameters were silently assumed.');
passed(1) = true;
fprintf('  PASS  %s\n',testNames{1});

%% Constant-coefficient SSPP analytic dispersion and exact monodromy
staticCfg = struct('topology','sspp','allowAssumptions',true, ...
    'a',5e-3,'Ls',45e-9,'mutualS',3e-9,'C0',20e-12, ...
    'deltaC',0,'Rs',0,'Gp',0,'fmHz',600e6);
staticModel = tl_build_model(staticCfg);
kStatic = [0.13 0.37 0.71]*pi/staticModel.cell.a;
Lk = staticModel.series.Ls+2*staticModel.series.mutualS* ...
    cos(kStatic*staticModel.cell.a);
omegaAnalytic = 2*abs(sin(kStatic*staticModel.cell.a/2))./ ...
    sqrt(Lk*staticModel.shunt.C0);
tmmStatic = tl_tmm_bands(staticModel,kStatic,struct('temporalSlices',32));
expected = [-omegaAnalytic;omegaAnalytic];
expected = fold_frequency(expected,staticModel.modulation.OmegaRadPerSec);
staticError = spectral_set_error(tmmStatic.omegaFolded,expected, ...
    staticModel.modulation.OmegaRadPerSec);
assert(staticError < 1e-11,'Static TMM error %.3g exceeds 1e-11.',staticError);
passed(2) = true;
fprintf('  PASS  %s (error %.3g)\n',testNames{2},staticError);

%% Modulated PWE/TMM cross-check
dynamicModel = tl_build_model(struct('topology','sspp', ...
    'allowAssumptions',true));
kDynamic = [0.17 0.39 0.63 0.84]*pi/dynamicModel.cell.a;
fourier = tl_pwe_fourier(dynamicModel,struct('Mtime',5,'Nt',256));
pwe = tl_pwe_bands(fourier,dynamicModel,kDynamic,struct());
tmm = tl_tmm_bands(dynamicModel,kDynamic,struct('temporalSlices',1024));
pweTmmError = spectral_set_error(pwe.omegaSelected,tmm.omegaFolded, ...
    dynamicModel.modulation.OmegaRadPerSec);
assert(pweTmmError < 0.01, ...
    'PWE/TMM circular spectral error %.3g exceeds 0.01.',pweTmmError);
passed(3) = true;
fprintf('  PASS  %s (error/Omega %.3g)\n',testNames{3},pweTmmError);

%% Passive loss sign under exp(i*k*x-i*omega*t)
lossModel = tl_build_model(struct('topology','sspp', ...
    'allowAssumptions',true,'deltaC',0,'Rs',0.4,'Gp',2e-4));
lossBands = tl_tmm_bands(lossModel,kDynamic,struct('temporalSlices',64));
passiveMaximum = max(imag(lossBands.omegaFolded),[],'all')/ ...
    lossModel.derived.omegaMaximumEstimate;
assert(passiveMaximum <= 1e-11, ...
    'A passive model has positive Im(omega): %.3g.',passiveMaximum);
passed(4) = true;
fprintf('  PASS  %s\n',testNames{4});

%% Finite chain energy and time-interface charge continuity
finiteModel = tl_build_model(struct('topology','sspp', ...
    'allowAssumptions',true,'cellCount',16,'deltaC',0));
nodeNumber = (0:finiteModel.finite.cellCount-1).';
initialVoltage = exp(1i*2*pi*nodeNumber/finiteModel.finite.cellCount);
fdtdCfg = struct('dt',0.05/finiteModel.derived.omegaMaximumEstimate, ...
    'nSteps',500,'recordEvery',5,'boundaryType','periodic', ...
    'modulationEnabled',false,'V0',initialVoltage);
field = tl_fdtd1d(finiteModel,fdtdCfg);
energyDrift = (max(field.energy.total)-min(field.energy.total))/ ...
    mean(field.energy.total);
assert(energyDrift < 1e-3,'FDTD energy drift %.3g exceeds 0.1%%.',energyDrift);

interfaceModel = tl_build_model(struct('topology','sspp', ...
    'allowAssumptions',true,'cellCount',8));
stepsPerPeriod = 64;
interfaceCfg = struct('dt',interfaceModel.modulation.period/stepsPerPeriod, ...
    'nSteps',stepsPerPeriod+2,'recordEvery',1,'boundaryType','periodic', ...
    'modulationEnabled',true, ...
    'modulationStart',interfaceModel.modulation.period,'V0',ones(8,1));
interfaceField = tl_fdtd1d(interfaceModel,interfaceCfg);
qBefore = interfaceField.node.Q(stepsPerPeriod,:);
qAt = interfaceField.node.Q(stepsPerPeriod+1,:);
chargeJump = norm(double(qAt-qBefore))/max(norm(double(qBefore)),realmin);
assert(chargeJump < 1e-12, ...
    'Charge changed across a uniform capacitance time interface.');
passed(5) = true;
fprintf('  PASS  %s (energy drift %.3g)\n',testNames{5},energyDrift);

%% CROW finite-Cblock finite-chain path
crowFinite = tl_build_model(struct('topology','crow', ...
    'allowAssumptions',true,'cellCount',8));
crowCfg = struct('dt',0.04/crowFinite.derived.omegaMaximumEstimate, ...
    'nSteps',40,'recordEvery',4,'boundaryType','matched', ...
    'modulationEnabled',false,'V0',[1;zeros(7,1)]);
crowField = tl_fdtd1d(crowFinite,crowCfg);
assert(isfield(crowField.resonator,'Qblock') && ...
    all(isfinite(crowField.resonator.Qblock),'all'), ...
    'CROW finite-Cblock state was not returned correctly.');
passed(6) = true;
fprintf('  PASS  %s\n',testNames{6});

%% Synthetic x-t transform at exact signed bins
a = 4e-3;
T = 1e-6;
nTime = 128;
nCell = 16;
dt = T/16;
time = (0:nTime-1).'*dt;
cells = 0:nCell-1;
x = cells*a;
Omega = 2*pi/T;
kExpected = 3*2*pi/(nCell*a);
omegaExpected = 0.25*Omega;
data = exp(-1i*omegaExpected*time)*exp(1i*kExpected*x);
observation = struct('data',data,'x',x,'t',time, ...
    'cellIndex',cells,'channelId',ones(1,nCell), ...
    'channelOffset',zeros(1,nCell),'observableName','synthetic I');
fft = tl_xt_fft_bands(observation,observation,struct( ...
    'cellPeriod',a,'temporalPeriod',T,'zeroPaddingTime',2, ...
    'zeroPaddingSpace',2,'referenceThresholdDb',80,'ridgeCount',1));
[~,peakLinear] = max(fft.foldedPower,[],'all','linear');
[peakOmegaIndex,peakKIndex] = ind2sub(size(fft.foldedPower),peakLinear);
assert(abs(fft.k(peakKIndex)-kExpected) < 1e-12*abs(kExpected));
assert(abs(fft.foldedOmega(peakOmegaIndex)-omegaExpected) < ...
    1e-12*abs(omegaExpected));
passed(7) = true;
fprintf('  PASS  %s\n',testNames{7});

%% Synthetic multi-probe time FFT with odd sampling and explicit folding
probeT = 1e-6;
probeOmega = 2*pi/probeT;
probeDt = probeT/5;
probePeriods = 7;
probeTime = (0:probePeriods*5).'*probeDt;
probeK = [0.2 0.6];
probeTone = (1+2/probePeriods)*probeOmega;
probeData = complex(zeros(numel(probeTime),2,numel(probeK)));
for index = 1:numel(probeK)
    carrier = exp(-1i*probeTone*probeTime);
    probeData(:,1,index) = carrier;
    probeData(:,2,index) = (0.4+0.2i)*carrier;
end
probeFft = tl_fdtd_fft_bands(probeData,probeTime,probeK,struct( ...
    'temporalPeriod',probeT,'analysisTimeRange',[0 probePeriods*probeT], ...
    'zeroPaddingFactor',1,'returnProbePower',true));
[~,probePeak] = max(probeFft.foldedPower,[],1);
assert(all(abs(probeFft.omegaOverOmega(probePeak)-2/probePeriods) < 1e-12));
[~,probeRawPeak] = max(probeFft.rawPower,[],1);
assert(all(abs(probeFft.rawOmegaOverOmega(probeRawPeak)- ...
    (1+2/probePeriods)) < 1e-12), ...
    'The full physical-frequency spectrum was not retained correctly.');
assert(probeFft.probeCount == 2 && ...
    size(probeFft.rawProbePower,2) == 2, ...
    'Fixed-probe power was not retained per probe.');
passed(8) = true;
fprintf('  PASS  %s\n',testNames{8});

%% Finite-chain Gaussian packets use complete-space FDTD probe data
packetModel = tl_build_model(struct('topology','sspp', ...
    'allowAssumptions',true,'cellCount',128,'deltaC',0, ...
    'fmHz',600e6));
packetT = packetModel.modulation.period;
packetK = [0.25 0.5 0.75]*pi/packetModel.cell.a;
packetScanCfg = struct('kScan',packetK,'dt',packetT/64, ...
    'nSteps',8*64,'recordEvery',1,'pulseIntensityFwhmCells',14, ...
    'voltageAmplitude',1,'centerNode',64, ...
    'probeNodeIndices',[59 62 64 66 69], ...
    'boundaryType','open','modulationEnabled',false, ...
    'initialTailTolerance',1e-4,'requireNoBoundaryArrival',true, ...
    'progressEvery',10);
packetScan = tl_fdtd_gaussian_k_scan(packetModel,packetScanCfg);
packetFft = tl_fdtd_fft_bands(packetScan.probeSignals, ...
    packetScan.time,packetK,struct('temporalPeriod',packetT, ...
    'analysisTimeRange',[0 8*packetT],'zeroPaddingFactor',2));
packetTmm = tl_tmm_bands(packetModel,packetK,struct('temporalSlices',32));
[~,packetPeak] = max(packetFft.foldedPower,[],1);
packetError = 0;
for index = 1:numel(packetK)
    distance = circular_distance(packetFft.foldedOmega(packetPeak(index)), ...
        real(packetTmm.omegaFolded(:,index)), ...
        packetModel.modulation.OmegaRadPerSec);
    packetError = max(packetError,min(distance)/ ...
        packetModel.modulation.OmegaRadPerSec);
end
assert(packetError <= 2*packetFft.nativeOmegaResolutionNormalized, ...
    'Finite-chain Gaussian FFT ridge exceeds two native frequency bins.');
assert(packetScan.initialEdgeAmplitudeBound < ...
    packetScan.initialTailTolerance, ...
    'Finite-chain Gaussian source does not satisfy its tail audit.');
passed(9) = true;
fprintf('  PASS  %s (error/Omega %.3g)\n',testNames{9},packetError);

%% Stroboscopic bulk FFT peak against a TMM eigenfrequency set
kBloch = [0.23 0.51 0.79]*pi/dynamicModel.cell.a;
bloch = tl_bloch_fft_bands(dynamicModel,kBloch,struct( ...
    'stepsPerPeriod',256,'temporalCellCount',64, ...
    'zeroPaddingFactor',4));
tmmBloch = tl_tmm_bands(dynamicModel,kBloch, ...
    struct('temporalSlices',1024));
blochError = 0;
for index = 1:numel(kBloch)
    distance = circular_distance(bloch.ridgeOmega(index), ...
        real(tmmBloch.omegaFolded(:,index)), ...
        dynamicModel.modulation.OmegaRadPerSec);
    blochError = max(blochError,min(distance)/ ...
        dynamicModel.modulation.OmegaRadPerSec);
end
assert(blochError <= 2*bloch.nativeOmegaResolution/ ...
    dynamicModel.modulation.OmegaRadPerSec, ...
    'Bloch FFT ridge error exceeds two native frequency bins.');
passed(10) = true;
fprintf('  PASS  %s (error/Omega %.3g)\n',testNames{10},blochError);

%% Fit and selection interfaces on synthetic/target-only data
fitK = [0.3 0.6].'*pi/staticModel.cell.a;
fitFrequency = zeros(size(fitK));
for index = 1:numel(fitK)
    [Aconstant,AinverseC] = staticModel.functions.bulkMatrices(fitK(index));
    fitFrequency(index) = max(real(1i*eig(Aconstant+ ...
        AinverseC/staticModel.shunt.C0)))/(2*pi);
end
fitData = struct('staticDispersion',struct('kRadPerM',fitK, ...
    'frequencyHz',fitFrequency,'uncertaintyHz',1e4));
fitBase = staticCfg;
fitBase.Ls = 40e-9;
fitResult = tl_fit_parameters(fitBase,fitData,struct( ...
    'parameterNames',{{'Ls'}},'maxIterations',40,'maxEvaluations',100));
assert_relative(fitResult.parameterValues.Ls,staticCfg.Ls,1e-6, ...
    'synthetic Ls fit');
selection = tl_select_components(staticModel,struct(),struct());
assert(selection.monteCarlo.sampleCount >= 1000 && ...
    strcmp(selection.status,'simulation-seed component screening'));
passed(11) = true;
fprintf('  PASS  %s\n',testNames{11});

%% Code Analyzer
mFiles = dir('*.m');
analyzerCount = 0;
for index = 1:numel(mFiles)
    issues = checkcode(mFiles(index).name,'-id');
    if ~isempty(issues)
        fprintf('  Code Analyzer: %s has %d issue(s).\n', ...
            mFiles(index).name,numel(issues));
    end
    analyzerCount = analyzerCount+numel(issues);
end
assert(analyzerCount == 0,'Code Analyzer reported %d issue(s).',analyzerCount);
passed(12) = true;
fprintf('  PASS  %s\n',testNames{12});

summary = struct('matlabVersion',version,'moduleDirectory',moduleDirectory, ...
    'testNames',{testNames},'passed',passed,'allPassed',all(passed), ...
    'pweTmmErrorNormalized',pweTmmError, ...
    'fdtdEnergyDrift',energyDrift, ...
    'finitePacketErrorNormalized',packetError, ...
    'blochErrorNormalized',blochError);
fprintf('Transmission-line validation passed: %d/%d tests.\n', ...
    sum(passed),numel(passed));
clear cleanupDirectory;
end

% -------------------------------------------------------------------------
function assert_relative(actual,expected,tolerance,label)
errorValue = abs(actual-expected)/max(abs(expected),realmin);
assert(errorValue <= tolerance,'%s relative error %.3g exceeds %.3g.', ...
    label,errorValue,tolerance);
end

function errorValue = spectral_set_error(actual,expected,Omega)
assert(isequal(size(actual),size(expected)));
errorValue = 0;
for index = 1:size(actual,2)
    nMode = size(actual,1);
    permutations = perms(1:nMode);
    best = Inf;
    for candidate = 1:size(permutations,1)
        selected = expected(permutations(candidate,:),index);
        realError = circular_distance(real(actual(:,index)), ...
            real(selected),Omega);
        imagError = abs(imag(actual(:,index))-imag(selected));
        best = min(best,max(hypot(realError,imagError))/Omega);
    end
    errorValue = max(errorValue,best);
end
end

function distance = circular_distance(a,b,period)
distance = abs(mod(a-b+period/2,period)-period/2);
end

function folded = fold_frequency(omega,Omega)
folded = mod(real(omega)+Omega/2,Omega)-Omega/2+1i*imag(omega);
end
