% Generate data
clear;

noise = 10/100; % percent
npatients = 30; % number of patients
bias = 1.1; % slope bias

% Regions per patient
territories = {'LAD','LCx','RCA'};
nterritories = length(territories);

% Patient states during measurement
states = {'Rest','Stress'};
nstates = length(states);

% Real flow values
restFlow = 0.7*(1+0.2*randn(npatients,nterritories));
stressFlow = 2*(1+0.2*randn(npatients,nterritories));

% Baseline data with noise
data1 = cat(3,  restFlow.*(1+noise*randn(npatients,nterritories)), stressFlow.*(1+noise*randn(npatients,nterritories)));
% Follow-up data with noise and a bias
data2 = bias * cat(3,  restFlow.*(1+noise*randn(npatients,nterritories)), stressFlow.*(1+noise*randn(npatients,nterritories)));

% BA plot paramters
tit = 'Flow Repeatability'; % figure title
gnames = {territories, states}; % names of groups in data {dimension 1 and 2}
label = {'Baseline Flow','Follow-up Flow','mL/min'}; % Names of data sets
corrinfo = {'n','SSE','r2','eq'}; % stats to display of correlation scatter plot
BAinfo = {'RPC(%)'}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
if 1
	colors = 'br'; % colors for the data sets
else
	colors = [0 0 1;...
		      1 0 0];
end
symbols = ''; % symbols for the data sets (default)

% Generat figure with symbols
[cr, fig, statsStruct] = BlandAltman(data1, data2,label,tit,gnames,corrinfo,BAinfo,limits,colors,symbols);

% Generate figure with numbers of the data points (patients)
symbols = 'Num';
BlandAltman(data1, data2,label,tit,gnames,corrinfo,BAinfo,limits,colors,symbols)

% Display statistical results that were returned from analyses
disp('Statistical results:');
disp(statsStruct);
