% Correlation - draws a correlation graph for two datasets.
%
% Correlation(data1, data2)
% Correlation(data1, data2,label) - 
% Correlation(data1, data2,label,tit,gnames)
% Correlation(data1, data2,label,tit,gnames,corrinfo) - specifies what
% information to display on the correlation chart as a cell of string in 
% order of top to bottom. The following codes are available:
%  - 'eq'     slope and intercept equation
%  - 'int%'   intercept as % of mean values
%  - 'r'      pearson r-value
%  - 'r2'     pearson r-value squared
%  - 'rho'    Spearman rho value
%  - 'SSE'    sum of squared error
%  - 'n'      number of data points used
%  * if not specified or empty, default is: {'eq';'r2';'SSE';'n'}
% Correlation(fig, ...) - specify a figure handle in which to
% display figure in which the correlation will be displyed 
% Correlation(ah, ...) - specify an axes which will be replaced by the 
% correlation axes.
% cr = Correlation(...) - return the coefficient of reproducibility
% (1.96*sd)
% [cr fig] = Correlation(...) - also return the figure handles
% [cr fig sstruct] = Correlation(...) - also return the structure of
% statistics for the analysis
%
% See also: BlandAltman

function [fig sstruct] = Correlation(varargin)

if isscalar(varargin{1}) && isequal(size(varargin{1}),[1 1]) && ishandle(varargin{1})
	shift = 1;
	fig = varargin{1};
else
	shift = 0;
	fig = [];
end
data1 = varargin{shift+1};
data2 = varargin{shift+2};
if nargin>=shift+3
	label = varargin{shift+3};
else
	label = '';
end
if nargin>=shift+4
	tit = varargin{shift+4};
else
	tit = '';
end
if nargin>=shift+5
	gnames = varargin{shift+5};
else
	gnames = '';
end
if nargin>=shift+6 && ~isempty(varargin{shift+6})
	if ischar(varargin{shift+6})
		corrinfo = {varargin{shift+6}};
	else
		corrinfo = varargin{shift+6};
	end
else
	corrinfo = {'eq';'r2';'SSE';'n'};
end

if iscell(label)
	if length(label)==1
		labelx = label{1};
		labely = label{1};
		labelm = label{1};
		labeld = ['\Delta ' label{1}];
	elseif length(label)==2
		labelx = label{1};
		labely = label{2};
		labelm = ['Mean ' label{1} ' & ' label{2}];
		labeld = [label{2} ' - ' label{1}];
	else
		labelx = [label{1} ' (' label{3} ')'];
		labely = [label{2} ' (' label{3} ')'];
		labelm = ['Mean ' label{1} ' & ' label{2} ' (' label{3} ')'];
		labeld = [label{2} ' - ' label{1} ' (' label{3} ')'];
	end	
else
	labelx = label;
	labely = label;
	labelm = label;
	labeld = ['\Delta ' label];
end

colors = 'rbgmc';
symb = 'sodp^v';
markersize = 4;

s = size(data1);
if ~isequal(s,size(data2));
	error('data1 and data2 must have the same size');
end

switch length(s)
	case 1
		s = [s 1 1];
	case 2
		s = [s 1];
	case 3
	otherwise
		error('Data have too many dimension');
end
n = s(1); % number of elements in each group
groups = numel(data1)/n;
data1 = reshape(data1, [numel(data1),1]);
data2 = reshape(data2, [numel(data2),1]);
mask = isfinite(data1) & isnumeric(data1) & isfinite(data2) & isnumeric(data2);

if isempty(fig)
	fig = figure;
	set(fig,'units','centimeters','position',[5 5 12 10],'color','w');
	cah = axes;
elseif strcmpi(get(fig,'type'),'figure')
	cah = axes;
elseif strcmpi(get(fig,'type'),'axes')
	cah = fig;
	fig = get(cah,'parent');
else
	error('What in tarnations is the handle that was passed to Bland-Altman????')
end
set(cah,'tag','Correlation Plot');

%% Correlation
hold(cah,'on');
for i=1:groups
	if s(3)==1
		marker = [colors(i) symb(1)];
	else
		marker = [colors(floor((i-1)/s(2))+1) symb(rem(i-1,s(2))+1)];
	end
	ph=plot(cah,data1((i-1)*n+(1:n)),data2((i-1)*n+(1:n)),marker);
	set(ph,'markersize',markersize);
end
% Linear regression
[polyCoefs, S] = polyfit(data1(mask),data2(mask),1);
[r p] = corrcoef(data1(mask),data2(mask)); r=r(1,2); p=p(1,2);
rho = corr(data1(mask),data2(mask),'type','Spearman');
N = sum(mask);
SSE = sqrt(sum((polyval(polyCoefs,data1(mask))-data2(mask)).^2)/(N-2));
a = axis(cah);
plot(cah,a(1:2), polyval(polyCoefs,a(1:2)),'-k');
if 0 % Add 95% CI lines
	xfit = a(1):(a(2)-a(1))/100:a(2);
	[yfit, delta] = polyconf(polyCoefs,xfit,S);
	h = [plot(cah,xfit,yfit+delta);...
		plot(cah,xfit,yfit-delta)];
	set(h,'color',[0.6 0.6 0.6],'linestyle','-');
end
corrtext = {};
for i=1:length(corrinfo)
	switch lower(corrinfo{i})
		case 'eq', corrtext = [corrtext; ['y=' num2str(polyCoefs(1),3) 'x+' num2str(polyCoefs(2),3)]];
		case 'int%', corrtext = [corrtext; ['intercept=' num2str(polyCoefs(2)/mean(data1+data2)*2*100,3) '%']];
		case 'r2', corrtext = [corrtext; ['r^2=' num2str(r^2,4)]];
		case 'r', corrtext = [corrtext; ['r=' num2str(r,4)]];
		case 'p', corrtext = [corrtext; ['p=' num2str(p,4)]];
		case 'rho', corrtext = [corrtext; ['rho=' num2str(rho,4)]];
		case 'sse', corrtext = [corrtext; ['SSE=' num2str(SSE,2)]];
		case 'n', corrtext = [corrtext; ['n=' num2str(N)]];
	end
end
text(a(1)+0.01*(a(2)-a(1)),a(3)+0.9*(a(4)-a(3)),corrtext,'parent',cah);
xlabel(cah,labelx); ylabel(cah,labely);


if nargout>1
	sstruct = struct('N',N,...
		'r',r,...
		'r2',r^2,...
		'SSE',SSE,...
		'rho',rho,...
		'Slope',polyCoefs(1),...
		'Intercept',polyCoefs(2));
end