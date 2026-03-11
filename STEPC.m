%% ============================
%  STEP C 
%  ============================
clear; clc;

% ---- Upload candidate set
useAprime = exist('STEP_Aprime_outputs.mat','file') == 2;
if useAprime
    SA = load('STEP_Aprime_outputs.mat','XS_aprime','FS_aprime');
    XS = SA.XS_aprime;  
    FS = SA.FS_aprime;
else
    SB = load('STEP_B_outputs.mat','X_S','F_S');
    XS = SB.X_S;
    FS = SB.F_S;
end

% ---- upload ref
RR = load('REF_global.mat','ref_global');     
ref_global = RR.ref_global;
fprintf('[STEP C] Candidate set size: %d\n', size(FS,1));
fprintf('STEP C: using fixed REF = [%.6g %.6g %.6g]\n', ref_global);

% Check:
if any((ref_global - FS) <= 0, 'all')
    warning('C: alcuni punti toccano/superano il REF fisso (MIN). EHV è comunque calcolabile.');
end


% ---- Parameters
Ntheta   = 500;          % draw of p(theta)
alphaCVaR= [0.90 0.95];  % CVaR
rng(1);                  % riproducibile

assert(exist('sample_theta','file')==2, 'Manca sample_theta.m nel path.');
assert(exist('obj_fun3','file')==2, 'Manca obj_fun3.m nel path.');

% ---- Results preallocation
nSol = size(XS,1);
nObj = size(FS,2);
Frob = nan(nSol, nObj, Ntheta);  

% ---- Pre-sampling: 
Thetas = sample_theta(Ntheta);    


parfor j = 1:Ntheta
    th = Thetas(j);               
    Fj = nan(nSol, nObj);
    for i = 1:nSol
        Fj(i,:) = obj_fun3(XS(i,:), th);
    end
    Frob(:,:,j) = Fj;
end
fprintf('Valutazioni completate: %d sol x %d obj x %d draws\n', nSol, nObj, Ntheta);

% ---- Prior-predictive-based Metrics ----
% 1) P(non-dom)
P_nondom = zeros(nSol,1);
for j = 1:Ntheta
    Fj = Frob(:,:,j);
    maskND = pf_min(Fj);
    P_nondom = P_nondom + double(maskND);
end
P_nondom = P_nondom / Ntheta;

% 2) Regret 
regret = zeros(nSol,1);
for j = 1:Ntheta
    Fj = Frob(:,:,j);
    fmin = min(Fj,[],1); fmax = max(Fj,[],1);
    Fn = (Fj - fmin) ./ max(fmax - fmin, 1e-12);  % 0=best,1=worst
    s = mean(Fn,2);
    rj = max(s) - s;
    regret = regret + rj;
end
regret = regret / Ntheta;

% 3) CVaR_alpha 
CVaR = struct();
for a = 1:numel(alphaCVaR)
    alpha = alphaCVaR(a);
    cvarMat = nan(nSol, nObj);
    for m = 1:nObj
        Y = squeeze(Frob(:,m,:));          
        q = quantile(Y, alpha, 2);        
        for i = 1:nSol
            tail = Y(i, Y(i,:) >= q(i));   
            cvarMat(i,m) = mean(tail);
        end
    end
    CVaR.(sprintf('a%02d',round(alpha*100))) = cvarMat; % es.: a90, a95
end

% 4) Expected Hypervolume 
HV_draw = zeros(Ntheta,1);
parfor j = 1:Ntheta
    Fj = Frob(:,:,j);
    Fj = Fj(pf_min(Fj),:);                 
    HV_draw(j) = hv3d_min(Fj, ref_global);
end
EHV_mean = mean(HV_draw);
EHV_sd   = std(HV_draw);
fprintf('Expected HV across prior draws: mean = %.6f, sd = %.6f\n', EHV_mean, EHV_sd);

% ---- PND-based robustness filtering with threshold δ ----
delta = 0.5;                      % P(non-dom) ≥ 0.5
keep_delta = P_nondom >= delta;

% ---- Composit Ranking (↑P_nonDom, ↓regret, ↓CVaR(95%))
cvar95 = CVaR.a95;                         
zP  = zscore(P_nondom);                    
zRg = -zscore(regret);                    
zC  = -zscore(mean(cvar95,2));             
score = normalize(zP + zRg + zC, 'range'); 

% ---- Top-10 ----
k = min(10, nSol);
[score_sorted, idx_sorted] = sort(score, 'descend');
topk_idx = idx_sorted(1:k);


% ---- Figure 1: heatmap P(non-dom) ----
[~, ord] = sort(score, 'descend');
figure('Color','w'); 
imagesc(P_nondom(ord)');
colormap(parula); colorbar; caxis([0 1]);
set(gca,'YTick',[],'XTick',1:nSol,'XTickLabel',ord, 'XTickLabelRotation',90);
title('P(non-dom) (ord. per Score)'); xlabel('Solution index (sorted)'); ylabel('P(non-dom)');
exportgraphics(gca,'C_heatmap_PnonDom.png','Resolution',300);


%% ===== helper =====
function mask = pf_min(F)
    n=size(F,1); mask=true(n,1);
    for i=1:n
        if ~mask(i), continue; end
        Fi=F(i,:);
        for j=1:n
            if j==i || ~mask(j), continue; end
            Fj=F(j,:);
            if all(Fj<=Fi) && any(Fj<Fi)
                mask(i)=false; break;
            end
        end
    end
end

function S2=prune2D(S)
    m=true(size(S,1),1);
    for i=1:size(S,1)
        if ~m(i), continue; end
        for j=1:size(S,1)
            if j==i||~m(j), continue; end
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

