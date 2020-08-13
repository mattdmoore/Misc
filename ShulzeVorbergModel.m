
clear all,clc

load('ShulzeVorbModelDATA.mat') % folder: ece.kaya/MATLAB 

data = data - 100*1000;         % subtract the 100 ms before stim 
IRIs = diff(data,1,2);          % IRIs       
asyn = data - 450000;           % asynchrony

%% data sample

I = IRIs(1,34:end)./1000;                                           % IRIs in miliseconds 
I = I';
N = length(I);

% tapping outlier removal
for i = 1:N
    if I(i) < 450 * 50/100 || I(i) > 450 * 150/100
        I(i) = [];
    end
end

%% WK model 
% ------------------------------ 
% I(j) = C(j) - D(j-1) + D(j)                                   % eq. (1)
% ------------------------------ 
% GAM_1 = lag_1 / lag_0                                         % lag 1 serial correlation 
% ------------------------------
% lag_1 = E[ (I(j) - mean(I)) * (I(j-1) - mean(I)) ]            % eq. (2), lag 1 covariance 
% lag_1 = E[ ( (C(j) - D(j-1) + D(j)) - mean(I)) * ...          % eq. (2), substitute I(j) with eq (1)  
%               ((C(j-1) - D(j-1-1) + D(j-1)) - mean(I)) ]      % eq. (2) & variance formula = E[X(j) - mean(X)]
% lag_1 = -var(D)
% lag_0 = E[ (I(j) - mean(I)) * (I(j) - mean(I)) ]              % eq. (3), lag 0 covariance 
% lag_0 = var(I)                                                % eq. (3) = variance
% lag_0 = var(C) + 2*var(D)                                     % when var(I) is distributed 
% GAM_1 = -var(D) / var(C) + 2*var(D)
% GAM_1 = - 1 / (2 + ( var(C) / var(D) ))                       % eq. (4), lag 1 serial corr
% ------------------------------  
% var(C) = lag_0 + 2*lag_1                                      % eq. (5), by eq. (3) 
% var(D) = -lag_1                                               % eq. (2)


% Long way
% =========================================================================
% % lag_1 = E[ (I(j) - mean(I)) * (I(j-1) - mean(I)) ] % eq. (2), lag 1 covariance 
% k = 1;
% for j = 2:N
%     lag_1(k,1) = (I(j) - mean(I)) * (I(j-1) - mean(I));         
%     k = k + 1;
% end
% lag_1 = sum(lag_1) / N; % divide by N-1 bcs we lose 1 datapoint

% % lag_0 = E[ (I(j) - mean(I)) * (I(j) - mean(I)) ] % eq. (3), lag 0 covariance (variance) 
% k = 1;
% for j = 2:N
%     lag_0(k,1) = (I(j) - mean(I)) * (I(j) - mean(I));           
%     k = k + 1;
% end
% lag_0 = sum(lag_0) / N;                                           

I_0 = I(1:end-1);
I_1 = I(2:end);

% lag_0 = E[ (I(j) - mean(I)) * (I(j) - mean(I)) ] % eq. (3), lag 0 covariance (variance)

lag_0 = sum((I_0 - mean(I_0)).*(I_0 - mean(I_0)))./length(I_0);
var(I_0)


%%
clear all,clc
x = randi(100,15,1);
x = x - mean(x);
for i=1:length(x)
    y(i) = x(i)*x(i);
end

var(x)
sum(y)./length(y)

%%
% VAR_C = lag_0 + 2*lag_1;
% VAR_D = -lag_1;

% Short way
% =========================================================================
motorvar = -cov(I(1:end-1),I(2:end));
motorvar = motorvar(1,2);
clockvar = var(I) + 2* cov( I(1:end-1) , I(2:end) );
clockvar = clockvar(1,2);


%% Semjen et al 2000 

% BASE MODEL: In = Tn + Mn+1 - Mn
    % MODEL w/ err correction: In = (Tn - a*An - b*An-1) + Mn+1 - Mn 

% Present:
    % In = inter-response-interval (ITI)
    % An = asynchrony (asyn = tap onset - metronome onset) 

% To be inferred:
    % Tn = timekeeper / brain clock >> mean(T) & std(T)
    % M = motor delay >> mean(M) & std(M)
    % a = fixed proportion to subtract from the last sync. error
    % b = fixed proportion to subtract from the next-to-the-last sync. error

%% ACVF = Compute autocovariance. https://www.mathworks.com/matlabcentral/fileexchange/24066-autocov-m

% Remove mean from x:
x = x - mean(x);
% For faster running time, we pre-allocate the output array:
acv = zeros(maxlag+1,1);
% Compute autocovariance:
for h = 0 : maxlag
   % Take matrix product of row vector and column vector to obtain a
   % scalar.  The row vector contains the first m-h elements of x; the
   % column vector contains the last m-h elements of x.
   acv(h+1)= x(1:m-h)' * x(1+h:m);
end
acv = acv / m;

%% another eq from https://atmos.washington.edu/~breth/classes/AM582/lect/lect7-notes.pdf
% lag-p ACVF which measures how strongly a time series is related with itself p samples later or earlier


%% fit ?    
p = [2 5]; 
X = 0:20; 
Yfcn = @(p,X) 1./(1+exp(-p(1).*(X-p(2)))); 

Y = Yfcn(p,X) + 0.1*randn(size(X)); 
	
[par fit]=fminsearch(@(p) norm(1./(1+exp(-p(1).*(X-p(2)))) -Y), [1,1])
figure(1)
plot(X, Y, 'p')
hold on
plot(X, Yfcn(par,X), '-r')
hold off
grid

% find the parameters where observed values are closest to expected values
% (function) 