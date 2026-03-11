%+% =======================
%  STEP A — Nominal MOO 
%  =======================

clearvars; clc;

%% 1) Nominal problem
theta0  = default_theta();                 
fun_nom = @(x) obj_fun(x, theta0);        

% Bounds (10 decision variables in [0,1])
lb = zeros(1,10);
ub = ones(1,10);

%% 2) Option
popSize = 500;
nGen    = 200;
opts_ga = optimoptions('gamultiobj', ...
    'PopulationSize',        popSize, ...
    'MaxGenerations',        nGen, ...
    'FunctionTolerance',     1e-6, ...
    'CrossoverFraction',     0.9, ...
    'MutationFcn',           @mutationadaptfeasible, ...
    'UseVectorized',         false, ...
    'Display',               'off');

%% 3) Multiple run per seed 
seeds  = 1:10;
X_runs = cell(numel(seeds),1);  
F_runs = cell(numel(seeds),1);   

pf = @(F) is_non_dominated_min(F);  

for i = 1:numel(seeds)
    rng(seeds(i));
    [X,F] = gamultiobj(fun_nom, 10, [],[],[],[], lb, ub, [], opts_ga);

    % Filter
    isND      = pf(F);
    X_runs{i} = X(isND,:);
    F_runs{i} = F(isND,:);
end

nd_sizes = cellfun(@(A) size(A,1), F_runs);
fprintf('[STEP A] Per-seed ND size (mean ± sd): %.1f ± %.1f\n', mean(nd_sizes), std(nd_sizes));

%% 4) Global reference point for HV
allF       = vertcat(F_runs{:});
ref_global = max(allF,[],1) + 0.10*abs(max(allF,[],1));    

if any((ref_global - allF) <= 0, 'all')
    warning('Reference HV troppo stretto, aumento margine al 20%%.');
    ref_global = max(allF,[],1) + 0.20*abs(max(allF,[],1));
end

%% 5) Candidate set 
XS = vertcat(X_runs{:});
FS = vertcat(F_runs{:});

tol = 1e-6;
[FS, ia] = dedup_rows(FS, tol);
XS       = XS(ia,:);

isNDc = pf(FS);
FS    = FS(isNDc,:);
XS    = XS(isNDc,:);

fprintf('Candidate set size (after dedup & ND): %d\n', size(FS,1));

%% 6) Metrics

% HV 
haveHV = exist('hv3d_min','file')==2;
HV   = nan(numel(seeds),1);
spr  = nan(numel(seeds),1);

for i = 1:numel(seeds)
    Fi = F_runs{i};              

    % Hypervolume 3D (min)
    if haveHV
        HV(i) = hv3d_min(Fi, ref_global);
    end
    
    fmin = min(Fi,[],1); fmax = max(Fi,[],1);
    Fn   = (Fi - fmin) ./ max(fmax - fmin, 1e-12);
    D    = pdist2(Fn, Fn);
    D(1:size(D,1)+1:end) = inf;   
    dmin = min(D,[],2);
    spr(i) = std(dmin);   
end

if haveHV
    fprintf('HV per seed (mean ± sd)   = %.6f ± %.6f\n', mean(HV), std(HV));
end
fprintf('Spread per seed (mean ± sd)= %.6f ± %.6f\n', mean(spr), std(spr));

%% --- helpers ---
function isND = is_non_dominated_min(F)
    n = size(F,1); isND = true(n,1);
    for i = 1:n
        if ~isND(i), continue; end
        dom = all(bsxfun(@le, F, F(i,:)), 2) & any(bsxfun(@lt, F, F(i,:)), 2);
        dom(i) = false;
        if any(dom), isND(i) = false; end
    end
end

function [A, ia] = dedup_rows(A, tol)
    if nargin<2, tol=0; end
    [n, ~] = size(A); keep = true(n,1); ia = zeros(0,1);
    for i = 1:n
        if ~keep(i), continue; end
        ia(end+1,1) = i; %#ok<AGROW>
        d = max(abs(A - A(i,:)), [], 2);
        dup = d <= tol; dup(1:i) = false; keep(dup) = false;
    end
    A = A(ia,:);
end

% Save output STEP A
save STEP_A_outputs.mat X_runs F_runs XS FS ref_global lb ub theta0 seeds HV spr

%% 7) Visualizzation

% Boxplot HV per seed
figure('Color','w'); 
boxplot(HV, 'Labels', {sprintf('n=%d', numel(HV))});
ylabel('Hypervolume (min, ref\_global)');
title('HV per seed (boxplot)');
grid on; box on;
saveas(gcf, 'STEPA_HV_boxplot.png');
print(gcf,'-depsc','-painters','STEPA_HV_boxplot.eps');

