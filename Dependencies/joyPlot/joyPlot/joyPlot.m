function varargout = joyPlot(data,x,offset,varargin)
%JOYPLOT Plot your data in a ridgeline representation
%   JOYPLOT(DATA,X,OFFSET) The array DATA should have a size m by n, where 
%   n is the number ofdatasets to plot and m is the sampling number. X is 
%   a vector containing the x-coordinates. OFFSET is a scalar that 
%   determines the displacement between plots. 
%
%   JOYPLOT(DATA,X,OFFSET,OVERLAPMETHOD) If the optional parameter
%   OVERLAPMETHOD is set to 'variable', OFFSET is a the overlap between 
%   datasets in percent. By default, OVERLAPMETHOD is set to 'constant'.
%
%   JOYPLOT(DATA,X,OFFSET,OVERLAPMETHOD,REVERSE) If REVERSE is set to
%   'true', the first row of DATA will be plotted at the bottom and the
%   last at the bottom. 
%
%   JOYPLOT(___,Name,Value) Specifies patch and line properties using one
%   or more 'Name,Value' pairs. For a list see Properties below.
%
%   [hf,hl] = JOYPLOT(___) Returns the patch and line handles.
%
%   [hs,hf,hl] = JOYPLOT(___) If the stroke is requested, its handle can be
%   returned as well.
%
%   [hs,hf,hl,hvl] = JOYPLOT(___) If additional vertical lines are plotted,
%   their handles can be returned too.
%
%   The function accepts the following Name,Value pairs:
%
%   - FaceColor: works just like the color in the PATCH function. It
%   additionally accepts the input 'position', which colors the faces
%   according to the position in the y-axis.
%   - FaceAlpha: works just like 'FaceAlpha' in the function PATCH.
%   - LineColor: works just like the color of function PLOT, with the
%   addition that it can be set to 'none' to delete the line.
%   - LineWidth: works just like 'LineWidth' in the function PLOT.
%   - StrokeColor: adds a stroke to each dataset with the specified color.
%   Here, the color works just like in the funciton PATCH.
%   - StrokeWidth: adds a stroke to each dataset with the specified width.
%   Again, this works like in the function PATCH.
%   - VLines: adds vertical lines in each dataset corresponding to the
%   x-coordinates provided. This is useful to plot median, mean or mode
%   values for example.
%   - VLinesColor: works just like the color of function PLOT.
%   - VLinesWidth: works just like 'LineWidth' in the function PLOT.
%
%   EXAMPLE
%   The classic joy plot from the cover of Joy Division's Unknown Pleasures 
%   by Peter Saville.
%   % Save the data from the web
%   filename = websave('pulsar.csv',['https://gist.githubusercontent.com/borgar/',...
%     '31c1e476b8e92a11d7e9/raw/0fae97dab6830ecee185a63c1cee0008f6778ff6/',...
%     'pulsar.csv']);
%   % Import it into MATLAB(R)
%   data = readmatrix(filename);
%   x = linspace(0,93,size(data,2));
%   % Create the figure and show the magic
%   figure
%   joyPlot(data',x,4)
%   set(gcf,'position',[500,100,560,680])
%   set(gca,'Visible','off', 'box','off','XTick',[],'YTick',[])
%
%   See also PLOT, PATCH, FILL.
%
%   Author: Santiago Benito, Ruhr-Universität Bochum
%   Contact: santiago.benito@rub.de
%   Last modified: 23.04.2020
%   Version: v0.1

%% Parse the inputs
% Initialize the input parser
p = inputParser;
% Set defaults
dOverlapMethod = 'constant';
dReverse = false;
dFaceColor = 'w';
dFaceAlpha = 1;
dLineColor = 'k';
dLineWidth = 1.5;
dStrokeColor = 'w';
dStrokeWidth = 3;
dVLines = 0;
dVLinesColor = '--k';
dVLinesWidth = 1;
% Validation function
expectedMethods = {'variable','constant'};
validOverlapMethod = @(x) isnumeric(x) ||...
    (ischar(x) && any(validatestring(x,expectedMethods)));
% Add required
addRequired(p,'data',@isnumeric);
addRequired(p,'x',@isnumeric);
addRequired(p,'offset',@isnumeric);
% Add optional parameters
addOptional(p,'overlapMethod',dOverlapMethod,validOverlapMethod)
addOptional(p,'reverse',dReverse,@islogical)
% Add pair parameters
addParameter(p,'FaceColor',dFaceColor)
addParameter(p,'FaceAlpha',dFaceAlpha,@isnumeric)
addParameter(p,'LineColor',dLineColor)
addParameter(p,'LineWidth',dLineWidth,@isnumeric)
addParameter(p,'StrokeColor',dStrokeColor)
addParameter(p,'StrokeWidth',dStrokeWidth,@isnumeric)
addParameter(p,'VLines',dVLines,@isnumeric)
addParameter(p,'VLinesColor',dVLinesColor)
addParameter(p,'VLinesWidth',dVLinesWidth,@isnumeric)
% Parse
parse(p,data,x,offset,varargin{:})
overlapMethod = p.Results.overlapMethod;
reverse = p.Results.reverse;
FaceColor = p.Results.FaceColor;
FaceAlpha = p.Results.FaceAlpha;
LineColor = p.Results.LineColor;
LineWidth = p.Results.LineWidth;
StrokeColor = p.Results.StrokeColor;
StrokeWidth = p.Results.StrokeWidth;
VLines = p.Results.VLines;
VLinesColor = p.Results.VLinesColor;
VLinesWidth = p.Results.VLinesWidth;

%% Data preparation
% If the user wants to return the handles, so be it
nout = max(nargout,1) - 1;
cnt = 0;

% Format x as a column vector
x = x(:);

% Extract some important parameters
mini = min(x);
maxi = max(x);
[m,n] = size(data);

% Apply the reverse option if needed
if ~reverse
    data = fliplr(data);
end

% Get the position of each dataset
if strcmp(overlapMethod,'variable')
    ypos = cumsum(max(data,[],1))*(1-offset);
    ypos = [0, ypos]; ypos(end) = [];
elseif strcmp(overlapMethod,'constant')
    ypos = 0:offset:offset*(n-1);
end

% Prepare the data for patch plotting
X = repmat([mini;x;maxi],1,n);
Y = data + repmat(ypos,m,1);
Y = fliplr([min(Y,[],1);Y;min(Y,[],1)]);


% Check if a colormap is to be imployed
if strcmp(FaceColor,'position')
    FaceColor = ypos;   
end

% Check if the stroke needs to be used
stroke = false;
if sum(strcmp(varargin,'StrokeColor') + ...
        strcmp(varargin,'StrokeWidth')) > 0
    stroke = true;
end

%% Patch & Plot
% Get the current axes, if given
ax = gca;
% Store the number of existing children
numelAx = numel(ax.Children);
hold on

% If a stroke is required, do it
if stroke
    Xs = X(2:end,:); Xs(end,:) = NaN(1,n);
    Ys = Y(2:end,:); Ys(end,:) = NaN(1,n);
    if strcmp(StrokeColor,'interp')
        cs = fliplr([data;NaN(1,n)]);
    else
        cs = FaceColor;
    end
    hs = fill(Xs,Ys,cs,'FaceAlpha',0,'EdgeColor',StrokeColor,...
        'LineWidth',StrokeWidth);
    if nout>0, cnt = cnt + 1; varargout{cnt} = hs; end
end
% Plot the rest of the stuff
hf = fill(X,Y,FaceColor,'FaceAlpha',FaceAlpha,'EdgeColor','none');
if nout>0, cnt = cnt + 1; varargout{cnt} = hf; end
if ~strcmp(LineColor,'none')
    hl = plot(X(2:end-1,:),Y(2:end-1,:),LineColor,'LineWidth',LineWidth);
    if nout>0, cnt = cnt + 1; varargout{cnt} = hl; end
end

% Additional vertical lines
if VLines ~=0
    yv = diag(interp1(x,data,VLines));
    XV = repmat(VLines,2,1);
    YV = [ypos;ypos + yv(:,1)'];
    hvl = plot(fliplr(XV),fliplr(YV),VLinesColor,'LineWidth',VLinesWidth);
    if nout>0, cnt = cnt + 1; varargout{cnt} = hvl; end
end

hold off

% Sort the elements of the plot handling the case that the figure already
% had some stuff in it.
numelAxAfter = numel(ax.Children) - numelAx;
ids = repmat(1:n,numelAxAfter/n,1) +...
    repmat((0:numelAxAfter/(n+1))'*n,1,n);
ids = ids(:)';
if numelAx > 0, ids = [ids,numelAxAfter+1:numelAxAfter+numelAx]; end
ax.Children = ax.Children(ids);

%% Y-Axis labels
yticks(ypos + diff([ypos,ypos(end)+max(data(:,end))])/2)
ylbls = num2cell(1:n);
if ~reverse, ylbls = fliplr(ylbls); end
yticklabels(ylbls)
xlim([mini,maxi])
set(ax,'TickLength',[0 0])

