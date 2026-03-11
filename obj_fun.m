function F = obj_fun(x, theta)

%VConstraint:
rate=zeros(1,10);
for i=1:10
    if (theta.k22(i)*(1- 0.2*x(6)+0.2*theta.k21-0.1*x(8)))/theta.k25(i) > 1
        rate(i) = 1;
    else
        rate(i) = (theta.k22(i)*(1- 0.2*x(6)+0.2*theta.k21-0.1*x(8)))/theta.k25(i);
    end
end

% aa: agricultural areas[crop]
aa = theta.k29;   
a = sum(aa);

% irr_wu: irrigation water use
scal_irr_wu = (theta.k9*(1-0.15*x(5)-0.05*theta.k15-0.15*x(6))*(2-x(7))*theta.k6); 
irr_wu = scal_irr_wu*aa.*(1-x(3)*theta.k11+x(4)*theta.k12).*theta.k13.*(1-theta.k14);

% capil: capillarity rise - volume
capil = theta.k19 .* theta.k13 .* aa;  

% aa sum / Initial crop area sum
rapp = a/sum(theta.k24);

rapp2 = (rate .* aa)/a;

GWa = theta.k0;
gwa_t = theta.k0;
SE = theta.k23;
SQ = theta.k20;

theta_k7  = theta.k7_month;  
theta_k10 = theta.k10_month; 


T=360;


for t=1:T

    mese = mod(t - 1, 12) + 1;
    theta.k7  = theta_k7(mese);
    theta.k10 = theta_k10(mese);

    
%% GW Availability - OBJECTIVE 1
    
    if gwa_t <= theta.k1
        GWA_ir = x(1)*theta.k2*(1+(theta.C*theta.k3)/12);
        outflow = 0;
    else 
        GWA_ir = 0;
        outflow = gwa_t - theta.k1;
    end
    
    GWA_dr = (((1+x(2))*theta.k4*theta.k5*theta.k6*theta.k7)*(1-theta.k8) ...
          + theta.k10*sum(irr_wu)*(1-theta.k16) ...
          + theta.k17*(1-theta.k18)) ...
          + theta.k6*theta.k10*sum(capil) + outflow;

    gwa_t = GWA_ir - GWA_dr;
    GWa = gwa_t + GWa;
     
    gwa_t =GWa;
    
%% Soil Erotion 
    if SE>1
        SE = SE -(0.001 *(x(4)+x(8)+x(6))*(1/3));
    elseif SE>0
        SE = SE + 0.001*(rapp +theta.C + (1-x(6))+(1-x(8)))*(1/4);
    else
        SE = SE + 0.001*(rapp +theta.C + (1-x(6))+(1-x(8)))*(1/4)-(0.001 *(x(4)+x(8)+x(6))*(1/3));
    end
    
%% SW quality - OBJECTIVE 2
    NloadArea = (1-0.22*x(4)*t/360)*sum(rapp2);
    SOM = ((1-NloadArea) + (x(6)+x(4))*(t/360))*1/3;
    SQir = SQ*0.005 * (x(8)+SOM)/2;
    if SQ <1
        SQdr = SQ*0.001*(NloadArea + theta.C+(1-SOM) + SE/4)*rapp;
    else 
        SQdr = 0;
    end
    SQ = SQ +SQir -SQdr;
    
    
    %% Average agricultural sustainability (AAS) - OBJECTIVE 3
    nload = (theta.k22./theta.k25)*(1-0.2*x(6)-0.2*theta.k21-0.1*x(8));
    % Agricultural productivity:
    AP = ones(1,10) .*(theta.k26*(theta.k31 ./ theta.k31)*(SOM + x(5) + x(6)*t/360 +SQ)*(1/4));
    % Crop profitability:
    CP = theta.k27 .*(theta.k30 ./ theta.k30).*(1+ 0.2*theta.k28 - 0.4*(ones(1,10)-nload));   %vettore
    scal = 1+0.3*x(9)-0.12*(1-x(8))-0.5*(0.6*((1-theta.k16)*(1-GWa/theta.k1) + (1-x(7))-x(3))*(1/3) + 0.2*(1-x(10))+0.2*(1-x(6)))*(1/3);
    %Agricultural sustainability:
    AS = 0.5*(AP + CP*scal);
    % Average agricultural sustainability
    AAS = (1/10)*sum(AS);
end
f1 = -GWa;
f2 = -(SQ);
f3 = -(AAS);

%%
F = [f1,f2,f3];
end