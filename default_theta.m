function theta = default_theta()
% Nominal values

theta.k0  = 25;          % initial GW storage
theta.k1  = 37;          % GW storage threshold
theta.C   = 1;
theta.k2  = 1;
theta.k3  = -0.005;      % yearly change in groundwater recharge rate
theta.k4  = 2723.93;     % population
theta.k5  = 9;
theta.k6  = 1e-6;
theta.k8  = 0;
theta.k9  = 1.06;

theta.k11 = [0,0,0,0,0,0,0,0,0.176,0];
theta.k12 = [0,0,0,0,0,0,0,0,0.05,0];

theta.k13 = [4200,4200,0,4813,3300,4903,2650,1130,10462,3662];  % Avg crop water need

theta.k19 = [0,0,0.51,0.55,0.51,0.55,1,0.51,0,0.51];            % Capillary rise potential
theta.k14 = 0.88 .* theta.k19;

theta.k15 = 0;
theta.k16 = 0.8;
theta.k17 = 0.042;       % industrial water demand
theta.k18 = 0;

theta.k20 = 0.5;
theta.k21 = 0;

theta.k22 = [43,43,52,250,43,8,112,43,43,220];      % baseline fertilizer rate
theta.k23 = 0.5;

theta.k24 = [1000,1000,706.9,398.9,50.6,863,70.8,112.7,596,136.3];
theta.k25 = [45,45,101,201,110,28,84,89,45,168];

theta.k26 = 1/1.1;
theta.k27 = 0.5;
theta.k28 = 0.5;

theta.k29 = [1.7,12.7,706.9,398.9,50.6,863,70.8,112.7,596,136.3]; % initial crops area

% Stagional
theta.k7_month  = [1 1 1 1 1 1.6 1.6 1.6 1 1 1 1];   % fluctuating population coeff
theta.k10_month = [0 0 0 0 0.07 0.2 0.3 0.3 0.13 0 0 0];  % irrigation demand coeff

% k30 Target crop profitability
theta.k30 = [27470, 4812, 4000, 11018, 1291, 3000, 3000, 3000, 24931, 43000];
% k31 Average unit crop yield
theta.k31 = [0.44, 1.35, 0.22, 0.20, 2.67, 0.01, 0.01, 0.01, 0.65, 0.01];

         
end



