%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PC4 - Econometría Financiera
%%%% Modelo SVAR - Restricciones de Corto Plazo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Integrantes
% Luis Quispe
clear all;
clc;
rng(1234);   
%Basado en códigos elaborados por J. Martinez (PUCP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
clc;
rng(1234);     

%% Importing Data

datos = xlsread('preg_1.xlsx','Mensuales');

%select data
tt     = datos(:,1);   % Terms of trade 
i_ib     = datos(:,6);   % Monetary Policy Interest Rate
bvl     = datos(:,3);   % BVL Stock Index
tc     = datos(:,4);   %exchange rate (PEN to USD)

D = datetime(2012,1,1):calmonths(1):datetime(2022,04,30)

%plot the series
subplot (2,2,1), plot(D,tt), title('Términos de intercambio');
subplot (2,2,2), plot(D,i_ib), title('Tasa de interés i.b.');
subplot (2,2,3), plot(D,bvl), title('Índice de la BVL');
subplot (2,2,4), plot(D,tc), title('Tipo de cambio');

% Looking at the graphs it is possible that all series have unit root.
% Also, give a brief description of the series 

%% Testing in-level series
%Unit Root Test
%Using ADF test in a trend stationary model:

adftest(tt)
adftest(i_ib)
adftest(bvl)
adftest(tc)
% all tests do not reject the null hypothesis

%Using Phillips-Perron test:

pptest(tt)
pptest(i_ib)
pptest(bvl)
pptest(tc)
% all tests do not reject the null hypothesis, again!

%% Calculate returns for all series
[bvl_r] = tick2ret(bvl);
[i_pm_r] = tick2ret(i_ib);
[tt_r] = tick2ret(tt);
[tc_r] = tick2ret(tc);

% plot the series
subplot (2,2,1), plot(tt_r), title('Términos de intercambio en retornos');
subplot (2,2,2), plot(i_pm_r), title('Tasa de interés i.b. en retornos');
subplot (2,2,3), plot(bvl_r), title('Índice de la BVL en retornos');
subplot (2,2,4), plot(tc_r), title('Tipo de cambio en retornos');

%% Test returns

%Using ADF test in a trend stationary model:

adftest(tt_r)
adftest(i_pm_r)
adftest(bvl_r)
adftest(tc_r)
% all tests do reject the null hypothesis

%Using Phillips-Perron test:

pptest(tt_r)
pptest(i_pm_r)
pptest(bvl_r)
pptest(tc_r)
% again, all tests reject the null hypothesis. 
% Now we have stationary series, continue estimating

%% Reduced VAR
%Proposed order is: tt -> tc -> i -> bvl

y1      = [tt_r tc_r bvl_r i_pm_r]; 
names1  = [' T.T.                    ';
           ' F.X.                    ';
           ' BVL                     ';
           ' Interbank rate          '];
order   = [1 2 3 4]; 

%% Parameters

lag             = 2;                % Lags of VAR model
c               = 2;                % 0: no constant, 1: constant, 2:linear trend
T               = size(y1,1);       % Size of data
n               = size(y1,2);       % Number of variables
k               = size(y1,2)*lag+c; % Number of regressors
S               = 12;               % Numbers of periods for IRFs and FEVD

%% Order of variables in VAR model

y       = [];
names   = ['                         '];

for i = 1:n  
    y           = [y y1(:,order(i))];
    names(i,:)  = string(names1(order(i),:));
end

%% VAR Estimation

Y = y(lag+1:T,:); % Dropping first "lags" observation of y vector.
X = [];

switch(c)
    case 0        % No deterministic components.
        X = [];    
    case 1        % Include a constant.
        X = ones(T-lag,1);        
    case 2        % Include a linear trend.
        for t = 1:T-lag            
            X(t,:) = [1 t];            
        end
end

% Lets generate X matrix for lags of vector y.

for i = 1:lag
    x{i} = y(lag+1-i:T-i,:);
    X = [X x{i}];
end

B_ols = (X'*X)^(-1)*(X'*Y); % OLS estimator of the model.

e = Y-X*B_ols; % Reduced errors.

sigmae = e'*e/(T-lag-k); % Variance matrix of reduced errors.

%% Getting companion matrix F

AA = B_ols(c+1:end,:); % Coefficients of lags. Excluding coefficients of deterministic components.

F=[];

for i = 1:lag
    A{i}=AA(1+n*(i-1):n*i,:); % Generating A1,...,Ap matrices.
    F = [F A{i}];  % Generating first matrix row of companion matrix F.
end
 
F = [F;eye(n*(lag-1)),zeros(n*(lag-1),n)]; % Generating others matrix rows of companion matrix F.
J = [eye(n) zeros(n,n*(lag-1))];

%% Stability of the model

lambda = eig(F); % Calculating eigenvalues of companion matrix F
modulo = (real(lambda).^2 + imag(lambda).^2 ).^(1/2); % Calculating modulus of each eigenvalue.

disp('Modulus of eigenvalues');
disp(modulo);
if max(abs(modulo))<1
    disp('The model satisfy stationary conditions.');
elseif max(abs(modulo))>=1
    disp('The model does not satisfy stationary conditions.');
end

% Graph of unit circle and inverse roots

tempx = [-1:0.01:1];
tempy = (1-tempx.^2).^(1/2);
rootx = real(lambda)';
rooty = imag(lambda)';

figure('Name','Unit Circle')
h = scatter(rootx,rooty);
set(h,'MarkerEdgeColor','Blue')
set(h,'MarkerFaceColor','Blue')
hold on
h = plot(tempx,tempy,tempx,-tempy);
set(h(1),'Color','Black')
set(h(2),'Color','Black')
xlim([(min(tempx)-0.5) (max(tempx)+0.5)])
ylim([(min(tempy)-1.2) (max(tempy)+0.2)])
grid on
hold off
title('Inverse Roots of Characteristic Polynominal')

clear tempx tempy rootx rooty;
%The model satisfy stationary conditions.

%% Shock of Structural Form - Cholesky

units_structural_shock = 1; % 0: one unit shock, 1: one standard deviation shock

P = chol(sigmae)';
u = inv(P)*e';


switch(units_structural_shock)
    case 0
        shocks = eye(n).*(diag(P).^(-1)); % one unit shock
    case 1
        shocks = eye(n); % one std deviation shock
end

for i = 1:n
    
    name_sh_struc(i,:) = strcat(' Structural Shock ',' ',num2str(i));
    
end

IRF_s  = zeros(n,S+1,n);
MSE_s  = zeros(n,S+1);
Ome_s  = zeros(n,S+1,n);
FEVD_s = zeros(n,S+1,n);
omega  = zeros(n,n);

% Generating IRFs and FEVD

for j = 0:S
    
    Phi = J*(F'^j)*J';
    Theta = Phi*P*shocks;
    IRF_s(:,j+1,:) = Theta;
    Ome_s(:,j+1,:) = Theta.^2;
    
    if j == 0
        MSE_s(:,j+1) = diag(Phi*sigmae*Phi');
    elseif j>=1
        MSE_s(:,j+1) = MSE_s(:,j) + diag(Phi*sigmae*Phi');
    end
   
end

Ome_s = cumsum(Ome_s,2);
MSE_s2 = (MSE_s).^(1/2);

for s = 0:S
    
    FEVD_s(:,s+1,:) = Ome_s(:,s+1,:)./MSE_s(:,s+1)*100;
    
end


% Graph of IRF of Structural Shocks

zeroline = zeros(1,S+1);
irf_axis = 0:S;

switch(units_structural_shock)
    case 0
        figure('Name','IRF of Structural Shocks - One Unit Shock')
    case 1
        figure('Name','IRF of Structural Shocks - One Std. Dev. Shock')
end

cont = 1;

for var = 1:n
    
    for sh = 1:n
        
        subplot(n,n,cont)
        h = plot(irf_axis,zeroline,irf_axis,IRF_s(var,:,sh));
        set(h(1),'Color','Black','LineWidth',0.5)
        set(h(2),'Color','Black','LineWidth',1.5)
        title(strcat('Response of ',' ',names(var,:),' to ',name_sh_struc(sh,:)))
        axis tight;grid on;
        ylim([(min(IRF_s(var,:,sh))-0.01) max(IRF_s(var,:,sh)+.01)])
        cont = cont+1;
        
    end
    
end

% Graph of FEVD

colors = [1 1 1;
          0 0 1;
          1 1 0;
          1 0 1];

switch(units_structural_shock)
    case 0
        figure('Name','FEVD of Structural Shocks - One Unit Shock')
    case 1
        figure('Name','FEVD of Structural Shocks - One Std. Dev. Shock')
end

cont = 1;

for var = 1:n
    
    subplot(1,n,cont)
    indx = zeros(S+1,n);
    indx(:,:) = FEVD_s(var,:,:);
    h = bar(indx,'stacked');
    for i = 1:n
        set(h(i),'FaceColor',colors(i,:))
    end
    axis tight
    title(strcat('Decomposition of ',' ',names(var,:)))
    legend(h,name_sh_struc,'Location','southoutside')
    cont = cont+1;
    
end
