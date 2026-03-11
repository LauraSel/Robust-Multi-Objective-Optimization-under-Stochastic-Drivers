function Thetas = sample_theta(N)
base = default_theta();
Thetas = repmat(base, N, 1);

for j = 1:N
    th = base;

    % --- k3: Normal(0,0.005), trunc [-0.01, 0.01]
    th.k3 = truncnorm_rnd(0, 0.005, -0.01, 0.01);

    % --- k4: Weibull(shape=42.5507, scale=2758.44), trunc [2603.5, 2848]
    th.k4 = truncweibull_rnd(42.5507, 2758.44, 2603.5, 2848);

    % --- k7: monthly sequence
    th.k7_month = ones(1,12);
    th.k7_month(6) = tri_rnd(0.8, 1.6, 2.4);
    th.k7_month(7) = tri_rnd(0.8, 1.6, 2.4);
    th.k7_month(8) = tri_rnd(0.8, 1.6, 2.4);

    % --- k10: monthly sequence
    th.k10_month = zeros(1,12);
    th.k10_month(5) = tri_rnd(0.035, 0.07, 0.105);
    th.k10_month(6) = tri_rnd(0.10,  0.20, 0.30);
    th.k10_month(7) = tri_rnd(0.15,  0.30, 0.45);
    th.k10_month(8) = tri_rnd(0.15,  0.30, 0.45);
    th.k10_month(9) = tri_rnd(0.065, 0.13, 0.195);

    % --- k13: triangular (min, mode=peak, max)
    th.k13 = zeros(1,10);
    th.k13(1) = tri_rnd(2100,   4200,  6300);
    th.k13(2) = tri_rnd(2100,   4200,  6300);
    th.k13(4) = tri_rnd(2406.5, 4813,  7219.5);
    th.k13(5) = tri_rnd(1650,   3300,  4950);
    th.k13(6) = tri_rnd(2451.5, 4903,  7354.5);
    th.k13(7) = tri_rnd(1325,   2650,  3975);
    th.k13(8) = tri_rnd(565,    1130,  1695);
    th.k13(9) = tri_rnd(5231,   10462, 15693);
    th.k13(10) = tri_rnd(831,    3662,  5493);
    th.k13(3)= 0;

    % --- k17: triangular
    th.k17 = tri_rnd(0.021, 0.042, 0.063);

    % --- k19: triangular
    th.k19 = base.k19; 
    th.k19(3) = tri_rnd(0.255, 0.51, 0.765);
    th.k19(4) = tri_rnd(0.275, 0.55, 0.825);
    th.k19(5) = tri_rnd(0.255, 0.51, 0.765);
    th.k19(6) = tri_rnd(0.275, 0.55, 0.825);
    th.k19(7) = unifrnd(0.5, 1.0);
    th.k19(8) = tri_rnd(0.255, 0.51, 0.765);
    th.k19(10) = tri_rnd(0.255, 0.51, 0.765);
    

    th.k14 = 0.88 * th.k19;

    % --- k22: Uniform
    th.k22 = zeros(1,10);
    th.k22(1) = unifrnd(36, 50);
    th.k22(2) = unifrnd(36, 50);
    th.k22(3) = unifrnd(44, 60);
    th.k22(4) = unifrnd(213, 288);
    th.k22(5) = unifrnd(36, 50);
    th.k22(6) = unifrnd(6, 9);
    th.k22(7) = unifrnd(95, 129);
    th.k22(8) = unifrnd(36, 50);
    th.k22(9) = unifrnd(36, 50);
    th.k22(10)= unifrnd(187, 253);

    % --- k30: triangular
    th.k30 = zeros(1,10);
    th.k30(1)  = tri_rnd(13735,   27470,   41205);
    th.k30(2)  = tri_rnd(2406,    4812,    7218);
    th.k30(3)  = tri_rnd(2000,    4000,    6000);
    th.k30(4)  = tri_rnd(5509,    11018,   16527);
    th.k30(5)  = tri_rnd(645.5,   1291,    19365);   
    th.k30(6)  = tri_rnd(1500,    3000,    4500);
    th.k30(7)  = tri_rnd(1500,    3000,    4500);
    th.k30(8)  = tri_rnd(1500,    3000,    4500);
    th.k30(9)  = tri_rnd(12465.5, 24931,   37396.5);
    th.k30(10) = tri_rnd(21500,   43000,   64500);

    % --- k31: normal truncated
    th.k31 = zeros(1,10);
    th.k31(1)  = truncnorm_rnd(0.44,   0.11,    0.22,   0.66);
    th.k31(2)  = truncnorm_rnd(1.35,   0.3375,  0.675,  2.025);
    th.k31(3)  = truncnorm_rnd(0.22,   0.055,   0.11,   0.33);
    th.k31(4)  = truncnorm_rnd(0.20,   0.050,   0.10,   0.30);
    th.k31(5)  = truncnorm_rnd(2.67,   0.6675,  1.335,  4.005);
    th.k31(6)  = truncnorm_rnd(0.01,   0.0025,  0.005,  0.015);
    th.k31(7)  = truncnorm_rnd(0.01,   0.0025,  0.005,  0.015);
    th.k31(8)  = truncnorm_rnd(0.01,   0.0025,  0.005,  0.015);
    th.k31(9)  = truncnorm_rnd(0.65,   0.1625,  0.325,  0.975);
    th.k31(10) = truncnorm_rnd(0.01,   0.0025,  0.005,  0.015);

    Thetas(j) = th;
end
end

% ===== helper: triangular =====
function x = tri_rnd(a, m, b)
% a<=m<=b
u = rand();
if u < (m - a) / (b - a)
    x = a + sqrt(u*(b-a)*(m-a));
else
    x = b - sqrt((1-u)*(b-a)*(b-m));
end
end

% ===== helper: truncated normal =====
function x = truncnorm_rnd(mu, sigma, a, b)
max_tries = 1e5;
for i=1:max_tries
    z = mu + sigma*randn();
    if z >= a && z <= b
        x = z; return;
    end
end
x = min(max(z, a), b);
end

% ===== helper: truncated Weibull =====
function x = truncweibull_rnd(k_shape, lambda_scale, a, b)
max_tries = 1e5;
for i=1:max_tries
    z = wblrnd(lambda_scale, k_shape);
    if z >= a && z <= b
        x = z; return;
    end
end
x = min(max(z, a), b);
end
