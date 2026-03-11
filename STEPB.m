%% STEP B 
clear; clc;

S = load('STEP_A_outputs.mat','F_runs','X_runs','FS','XS','seeds');

R = load('REF_global.mat','ref_global');
ref_global = R.ref_global;                
F_runs = S.F_runs;                        
X_runs = S.X_runs;                        
FS     = S.FS;                             
XS     = S.XS;                             
seeds  = S.seeds;


% 1) Statistics
nd_sizes = cellfun(@(A) size(A,1), F_runs);
fprintf('Per-seed ND size (mean ± sd): %.2f ± %.2f\n', mean(nd_sizes), std(nd_sizes));

% 2) Union and dedup 
tol    = 1e-6;
Fall   = vertcat(F_runs{:});   
Xall   = vertcat(X_runs{:});

[FSu, XSu] = dedup_FX(Fall, Xall, tol);

% 3) Global ND 
isNDg = pf(FSu);
F_S   = FSu(isNDg,:);
X_S   = XSu(isNDg,:);
fprintf('Candidate set size (after dedup & ND): %d\n', size(F_S,1));

% 4) Reference check
Fall = vertcat(F_runs{:});
if any((ref_global - Fall) <= 0, 'all')
    warning('Reference fisso risulta stretto per alcuni punti (non lo modifico).');
end

% 5) HV per seed
HV_seed = nan(numel(seeds),1);
for i=1:numel(seeds)
    Fi = F_runs{i};
    Fi = Fi(pf(Fi),:);
    HV_seed(i) = hv3d_min(Fi, ref_global);   % MIN
end
fprintf('HV per seed (mean ± sd) = %.6f ± %.6f\n', mean(HV_seed), std(HV_seed));

% 6) Save output STEP B
save('STEP_B_outputs.mat','F_S','X_S','ref_global','HV_seed','nd_sizes','F_runs','X_runs');


%%
function mask = pf(F)
    n=size(F,1); mask=true(n,1);
    for i=1:n
        if ~mask(i), continue; end
        Fi=F(i,:);
        for j=1:n
            if j==i||~mask(j), continue; end
            Fj=F(j,:);
            if all(Fj<=Fi) && any(Fj<Fi)
                mask(i)=false; break;
            end
        end
    end
end

function [Fo,Xo] = dedup_FX(F,X,tol)
    if nargin<3, tol=1e-6; end
    sc = 1/tol;
    key = round([F X]*sc)/sc;
    [~,ia] = unique(key,'rows','stable');
    Fo = F(ia,:); Xo = X(ia,:);
end



function S2=prune2D(S)
    m=true(size(S,1),1);
    for i=1:size(S,1)
        if ~m(i), continue; end
        for j=1:size(S,1)
            if j==i || ~m(j), continue; end
            if all(S(j,:)>=S(i,:)) && any(S(j,:)>S(i,:)), m(i)=false; break; end
        end
    end
    S2=S(m,:);
end
function A=area2D_max(S)
    if isempty(S), A=0; return; end
    [~,idx]=sort(S(:,1),'descend'); S=S(idx,:);
    A=0; g2p=0; h=0;
    for i=1:size(S,1)
        g2=S(i,1); g3=S(i,2);
        w=max(g2-g2p,0); h=max(h,g3);
        A=A + w*h; g2p=g2;
    end
end


%% ===== EAF 3D (MIN) — L50 / L90 

if exist('F_runs','var')~=1
    tmp = load('STEP_A_outputs.mat','F_runs');
    F_runs = tmp.F_runs;
end

objNames = {'f1','f2','f3'};      
levels   = [0.5, 0.9];
outdir   = 'EAF3D_figs'; if ~exist(outdir,'dir'); mkdir(outdir); end

[L50, L90, info] = eaf3d_L_levels_array(F_runs, levels);

fprintf('[EAF3D] nRuns=%d | t50=%d, |L50|=%d | t90=%d, |L90|=%d\n', ...
    info.nRuns, info.tIdx(1), size(L50,1), info.tIdx(2), size(L90,1));

% Plot scatter3
figure('Color','w'); hold on; grid on; box on;
if ~isempty(L90)
    scatter3(L90(:,1), L90(:,2), L90(:,3), 26, 'filled', 'MarkerFaceAlpha', 0.55);
end
if ~isempty(L50)
    scatter3(L50(:,1), L50(:,2), L50(:,3), 36, 'filled', 'MarkerFaceAlpha', 0.90);
end
xlabel(sprintf('%s (min)', objNames{1}));
ylabel(sprintf('%s (min)', objNames{2}));
zlabel(sprintf('%s (min)', objNames{3}));
title(sprintf('EAF 3D — L50 (t=%d) e L90 (t=%d)', info.tIdx(1), info.tIdx(2)));
legend({'L90','L50'}, 'Location','best');

fn = fullfile(outdir, 'EAF3D_L50_L90.png');
exportgraphics(gca, fn, 'Resolution', 300);
fprintf('Saved: %s\n', fn);



%% ===== helper =====
function [L50, L90, info] = eaf3d_L_levels_array(F_runs, levels)

    if nargin<2 || isempty(levels), levels = [0.5 0.9]; end
    levels = levels(:)';
    nRuns  = numel(F_runs);
    tIdx   = max(1, min(nRuns, ceil(levels * nRuns)));

    for r=1:nRuns
        F = F_runs{r};
        if isempty(F), continue; end
        F_runs{r} = F(pf_min(F),:);
    end

    Zall = vertcat(F_runs{:});
    if isempty(Zall)
        L50 = []; L90 = [];
        info = struct('nRuns',nRuns,'levels',levels,'tIdx',tIdx,'nZ',0);
        return;
    end
    Zvals = unique(Zall(:,3), 'sorted');
   
    Lcell = cell(1, numel(tIdx));
    for k=1:numel(tIdx)
        t = tIdx(k);
        Lcell{k} = compute_Lt_3D_via_slices(F_runs, t, Zvals);
        Lcell{k} = Lcell{k}(pf_min(Lcell{k}),:); % ND finale 3D
    end

    L50 = Lcell{1};
    L90 = Lcell{2};

    info = struct('nRuns',nRuns,'levels',levels,'tIdx',tIdx,'nZ',numel(Zvals));
end


function Lt3 = compute_Lt_3D_via_slices(F_runs, t, Zvals)

    nRuns = numel(F_runs);
    Xall = vertcat(F_runs{:});
    xVals = unique(Xall(:,1), 'sorted');
    Kx = numel(xVals);

    R = cell(nRuns,1);
    for r=1:nRuns
        F = F_runs{r};
        if isempty(F)
            R{r} = zeros(0,3);
        else
            [~,ix] = sort(F(:,3),'ascend');
            R{r} = F(ix,:);
        end
    end

    Lt3_acc = zeros(0,3);

    for iz = 1:numel(Zvals)
        Z = Zvals(iz);

        yminMat = inf(nRuns, Kx);

        for r = 1:nRuns
            F = R{r};
            if isempty(F), continue; end

            last = find(F(:,3) <= Z, 1, 'last');
            if isempty(last), continue; end
            Fz = F(1:last, 1:2); % (x,y)
            Fz = Fz(pf_min_2d(Fz),:);

            [xs,ord] = sort(Fz(:,1),'ascend');
            ys = Fz(ord,2);
            ypref = cummin(ys); 
            idx = arrayfun(@(x) find(xs<=x,1,'last'), xVals, 'UniformOutput', false);
            idx = cellfun(@(u) iff(isempty(u),0,u), idx);

            yv = inf(1,Kx);
            ok = idx>0;
            yv(ok) = ypref(idx(ok));
            yminMat(r,:) = yv;
        end

       
        y_t = kth_smallest_per_column(yminMat, t);
        y_t = enforce_monotone_nonincreasing(y_t);
 
        y_t = y_t(:);          
        xv  = xVals(:);       
        
        chg  = [true; abs(diff(y_t))>0];
        keep = chg & isfinite(y_t);
        
        P2 = [xv(keep), y_t(keep)];

        if ~isempty(P2)
            P3 = [P2, Z*ones(size(P2,1),1)];
            Lt3_acc = [Lt3_acc; P3]; 
        end
    end

    Lt3 = Lt3_acc;
    if ~isempty(Lt3)
        Lt3 = Lt3(pf_min(Lt3),:);
    end
end


function yk = kth_smallest_per_column(M, k)
    [n, K] = size(M);
    if k>n, yk = inf(1,K); return; end
    Ms = sort(M, 1, 'ascend');
    yk = Ms(k,:);
end


function y = enforce_monotone_nonincreasing(y)
    y = y(:)';
    for i=2:numel(y)
        if y(i) > y(i-1)
            y(i) = y(i-1);
        end
    end
end


function m = pf_min(F)
    n=size(F,1); m=true(n,1);
    for i=1:n
        if ~m(i), continue; end
        Fi=F(i,:);
        for j=1:n
            if j==i||~m(j), continue; end
            Fj=F(j,:);
            if all(Fj<=Fi) && any(Fj<Fi)
                m(i)=false; break;
            end
        end
    end
end

function m = pf_min_2d(F2)
    n=size(F2,1); m=true(n,1);
    for i=1:n
        if ~m(i), continue; end
        Fi=F2(i,:);
        for j=1:n
            if j==i||~m(j), continue; end
            Fj=F2(j,:);
            if all(Fj<=Fi) && any(Fj<Fi)
                m(i)=false; break;
            end
        end
    end
end

function out = iff(cond, a, b)
    if cond, out = a; else, out = b; end
end



