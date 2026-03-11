% MAKE_REF.m 
clear; clc;

Fall = [];
src  = '';



S = load('STEP_A_outputs.mat');
if isfield(S,'F_runs') && ~isempty(S.F_runs)
    Fall = vertcat(S.F_runs{:});   
    src  = 'STEP_A_outputs.mat : vertcat(F_runs{:})';
elseif isfield(S,'FS') && ~isempty(S.FS)
    Fall = S.FS;                   
    src  = 'STEP_A_outputs.mat : FS';
end


% ---- Fixed reference ----
margin = 0.30;                         
Fmax   = max(Fall,[],1);
ref_global = Fmax + margin*abs(Fmax);

viol = any((ref_global - Fall) <= 0, 'all');
fprintf('Source: %s\n', src);
fprintf('REF = [%.6g  %.6g  %.6g]\n', ref_global);


save('REF_global.mat','ref_global');
fprintf('Salvato: REF_global.mat\n');

% ===== helper =====
function m = pf_min(F)
    n=size(F,1); m=true(n,1);
    for i=1:n
        if ~m(i), continue; end
        Fi=F(i,:);
        for j=1:n
            if j==i || ~m(j), continue; end
            Fj=F(j,:);
            if all(Fj<=Fi) && any(Fj<Fi)
                m(i)=false; break;
            end
        end
    end
end

