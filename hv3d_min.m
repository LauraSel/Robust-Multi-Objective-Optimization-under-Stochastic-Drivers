function HV = hv3d_min(F, ref)
% HV of a NON-DOMINATED set F (Kx3) for MINIMIZATION, w.r.t. reference point ref (1x3).
% Assumes ref >= each coordinate of F (worse point). Returns union volume of
% axis-aligned boxes [f1,ref1]x[f2,ref2]x[f3,ref3].

if isempty(F)
    HV = 0; return;
end
if size(F,2) ~= 3
    error('hv3d_min: F must be Kx3');
end
if numel(ref) ~= 3
    error('hv3d_min: ref must be 1x3');
end


mask = all(F <= ref, 2);
F = F(mask,:);
if isempty(F)
    HV = 0; return;
end


F = sortrows(F, 1);    % [f1 f2 f3]
K = size(F,1);
HV = 0;

for i = 1:K
    f1_i   = F(i,1);
    if i < K
        f1_next = F(i+1,1);
    else
        f1_next = ref(1);
    end
    width = f1_next - f1_i;
    if width <= 0
        continue;
    end

    
    P = F(1:i, 2:3); % [f2 f3]

   
    Pnd = nondom2d_min(P);

   
    area2D = area_union_rect2d_min(Pnd, ref(2), ref(3));

    HV = HV + width * area2D;
end
end

% -------------------- helpers --------------------

function Pnd = nondom2d_min(P)
N = size(P,1);
keep = true(N,1);
for i = 1:N
    if ~keep(i), continue; end
    pi = P(i,:);
    % qualcuno domina i?
    dom = all(bsxfun(@le, P, pi), 2) & any(bsxfun(@lt, P, pi), 2);
    dom(i) = false;
    if any(dom)
        keep(i) = false;
    end
end
Pnd = P(keep,:);
end

function A = area_union_rect2d_min(P, ref2, ref3)
if isempty(P)
    A = 0; return;
end

P = sortrows(P, [1 2]);

prev_f3 = ref3;
A = 0;
for i = 1:size(P,1)
    f2 = min(P(i,1), ref2);
    f3 = min(P(i,2), ref3);
    if f3 < prev_f3
        A = A + max(ref2 - f2, 0) * max(prev_f3 - f3, 0);
        prev_f3 = f3;
    end
end
end
