function varargout = BeadGazer2(varargin)
% BEADGAZER2 MATLAB code for BeadGazer2.fig
%      BEADGAZER2, by itself, creates a new BEADGAZER2 or raises the existing
%      singleton*.
%
%      H = BEADGAZER2 returns the handle to a new BEADGAZER2 or the handle to
%      the existing singleton*.
%
%      BEADGAZER2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEADGAZER2.M with the given input arguments.
%
%      BEADGAZER2('Property','Value',...) creates a new BEADGAZER2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BeadGazer2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BeadGazer2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BeadGazer2

% Last Modified by GUIDE v2.5 02-Aug-2021 14:35:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BeadGazer2_OpeningFcn, ...
                   'gui_OutputFcn',  @BeadGazer2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before BeadGazer2 is made visible.
function BeadGazer2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BeadGazer2 (see VARARGIN)

global beadGazerOpt
evalin('base','global beaddata beadGazerOpt beadSelection');
if ~isfield(beadGazerOpt,'logZero')
    beadGazerOpt.logZero = 0.01;
end
if ~isfield(beadGazerOpt,'chosenPositions')
    beadGazerOpt.chosenPositions = [];
end
% Choose default command line output for BeadGazer2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BeadGazer2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BeadGazer2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function eSetRange_Callback(hObject, eventdata, handles)

function eSetsToPlot_Callback(hObject, eventdata, handles)
plotROIs(handles);

function cbUseValues_Callback(hObject, eventdata, handles)
plotROIs(handles);

function cbPlotLogX_Callback(hObject, eventdata, handles)
plotROIs(handles);

function pbReadFile_Callback(hObject, eventdata, handles)
global state beadsurv beaddata beadSurveyorActive beadGazerData roiImages trainedClassifier
% need to be careful not to overwrite live data if BeadSurveyor is running now
if ~isempty(beadSurveyorActive) && beadSurveyorActive
    error('BeadSurveyor appear to be active [global beadSurveyorActive is true] -- will not read file because it risks destroying live data');
    return;
end
if isempty(beadGazerData) || ~isfield(beadGazerData,'pathname')
    try 
        fname = fullfile(state.files.savePath,'*_beaddata*.mat');
    catch
        fname = '*_beaddata*.mat';
    end
else
    fname = fullfile(beadGazerData.pathname,beadGazerData.filename);
end
[filename,pathname] = uigetfile('*_beaddata*.mat','BeadGazer:Read BeadData File',fname);
if filename==0, return; end
load([pathname filename]);
beadGazerData.pathname = pathname;
beadGazerData.filename = filename;
handles.tFilename.String = [pathname filename];
try
    nsets = size(beaddata.roiData,2);
    handles.eSetsToPlot.String = ['1:' num2str(nsets)];
    handles.lbColToPlot.String = beaddata.valueNames;
    handles.lbColToPlot.Value  = 1;
end
showSetList(handles);
fname = fullfile(beadGazerData.pathname,'roiImages');
[fname,pname] = uigetfile('roiImages.mat','BeadGazer:Read roiImages (usually in \\Set)',fname);
if fname==0, return; end
load([pname fname]);
try
    load('gyBeadClassifier.mat');
    trainedClassifier = trainedModel;
catch
    return
end
if ~isempty(roiImages) && ~isnumeric(trainedClassifier)
    beadImageAnalysis(roiImages,trainedClassifier);
end

function pbReplot_Callback(hObject, eventdata, handles)
global beaddata beadGazerOpt
try
    data = beaddata.roiData;
catch
    error('No valid roiData'); return
end
try
    sets = eval(handles.eSetRange.String);
catch
    sets = 1:size(data,2);
end
nphots = mean(data(:,sets,3),2); % mean over sets
sw = isIntensity(handles); % get the switch value F=lifetime, T=intensity

if handles.cbLastMFirst.Value
    ltimeLo = getVals(sw,data,sets(1));   % was data(:,sets(1),1);    % LT(,1) or intensity(,2)
    ltimeHi = getVals(sw,data,sets(end)); % was data(:,sets(end),1);  % LT(,1) or intensity(,2)
else
    if handles.cbDTBasedOnVal.Value % new method
        % looking for biggest difference between min and max concs
        % get the values/concentrations used for each set
        
        % when plotting values, valIdx says which column to use
        valIdx = handles.lbColToPlot.Value;
        nomValues =  beaddata.setInfo(sets,valIdx)';
        % get the lowest and highest values
        minVal  = min(nomValues);
        maxVal  = max(nomValues);
        % get all set numbers matching each value
        minSets = sets(nomValues==minVal);
        maxSets = sets(nomValues==maxVal);
        % get the mean outcome for min and max values (covers repeated values)
        ltimeLo = getVals(sw,data,minSets,'mean'); % was mean(data(:,minSets,1),2);    % LT(,1) or intensity(,2)
        ltimeHi = getVals(sw,data,maxSets,'mean'); % was mean(data(:,maxSets,1),2);    % LT(,1) or intensity(,2)    
    else
        ltimeLo = getVals(sw,data,sets,'min'); % was min(data(:,sets,1),[],2); % min over sets % LT(,1) or intensity(,2)
        ltimeHi = getVals(sw,data,sets,'max'); % was max(data(:,sets,1),[],2); % max over sets % LT(,1) or intensity(,2)
    end
end
if numel(sets)==1
    deltas = ltimeLo;
else
    if sw
        deltas = log10(ltimeHi ./ ltimeLo);  % log10(F/F0)
        deltas = min(max(deltas,-1),1); % limit to +/- 1 (tenfold range)
    else
        deltas = ltimeHi-ltimeLo;     % LThi - LTlo
    end
end

% default correspondence:
xvals   = nphots;
xtext   = '# photons in ROI';
yvals   = deltas;
if sw
    if numel(sets)>1
        ytext = 'log10(F/F0)';
    else
        ytext = 'F';
    end
else
    if numel(sets)>1
        ytext = 'deltaTau';
    else
        ytext = 'tau';
    end
end
roinums = 1:numel(nphots);

% flexible interpretation
if handles.cbCriterion.Value || handles.cbAltX.Value || handles.cbAltY.Value 
    nsets  = size(beaddata.setInfo,1);
    ncols  = sum(cellfun(@(c) ~isempty(c), beaddata.valueNames));
    cnames = beaddata.valueNames(1:ncols);
    setstr = cell(1,ncols);
    % assign s1..sN as the lifetime for each set
    for k=1:nsets
        eval(['s' num2str(k) ' =  beaddata.roiData(:,k,1);']); 
        setstr{k} = '[';
        for jj=1:ncols
            setstr{k} = [setstr{k} cnames{jj} '=' num2str(beaddata.setInfo(k,jj)) ' '];
        end
        setstr{k}(end) = ']';
    end
    % assign the classifier values (if available) as cs
    if isfield(beaddata,'cs')
        cs = beaddata.cs(:);
    end
    % define i1..iN as intensity for the set
    for k=1:nsets
        eval(['i' num2str(k) ' =  beaddata.roiData(:,k,2);']);
    end
    % define n1..nN as nphot for the set
    for k=1:nsets
        eval(['n' num2str(k) ' =  beaddata.roiData(:,k,3);']);
    end
    % now interpret and apply the user-entered formulas
    if handles.cbAltX.Value
        try
            [~,xvals] = evalc(handles.eAltX.String);
            handles.eAltX.BackgroundColor = 'w';
        catch
            handles.eAltX.BackgroundColor = [1 .5 .5];
            return
        end
        if sw
            xtext = ['F: ' handles.eAltX.String];  % intensity (probably)
        else
            xtext = ['LT: ' handles.eAltX.String]; % lifetime
        end
        for k=1:nsets
            xtext = strrep(xtext,['s' num2str(k)],setstr{k});
            xtext = strrep(xtext,['i' num2str(k)],setstr{k});  % just replace both
        end
    end
    if handles.cbAltY.Value
        try
            [~,yvals] = evalc(handles.eAltY.String);
            handles.eAltY.BackgroundColor = 'w';
        catch
            handles.eAltY.BackgroundColor = [1 .5 .5];
            return
        end
        if sw
            ytext = ['F: ' handles.eAltY.String]; % intensity (probably)
        else
            ytext = ['LT: ' handles.eAltY.String]; 
        end
        for k=1:nsets
            ytext = strrep(ytext,['s' num2str(k)],setstr{k});
            ytext = strrep(ytext,['i' num2str(k)],setstr{k});  % just replace both
        end
    end
    if handles.cbCriterion.Value
        try
            [~,sel] = evalc(handles.eCriterion.String);
            handles.eCriterion.BackgroundColor = 'w';
        catch
            handles.eCriterion.BackgroundColor = [1 .5 .5];
            return
        end
        xvals = xvals(sel);
        yvals = yvals(sel);
        roinums = roinums(sel);
    end
end

% set the colors if there has been a classification
if isfield(beaddata,'cs')
    if ~exist('sel','var'), sel=1:numel(xvals); end
    cols = [0 0 1; 0 .6 0];
    cs = cols((1+beaddata.cs(sel)),:);
else
    cs = 'b';
end


% restrict to an area of the field of view (don't use with criterion!)
if handles.cbPosFilter.Value && ~isempty(beadGazerOpt.chosenPositions)
    xvals   = xvals(beadGazerOpt.chosenPositions,:,:);
    yvals   = yvals(beadGazerOpt.chosenPositions,:,:);
    roinums = roinums(beadGazerOpt.chosenPositions);
    cs = 'b';
end

% plot the x-histogram
axes(handles.axHistNPhot);
hist(xvals,50);
xlabel(xtext);
ylabel('# ROIs');
% plot the y-histogram
axes(handles.axHistDeltaLT);
hist(yvals,50);
set(gca,'View',[90 -90]);
xlabel(ytext);
ylabel('# ROIs');
% plot the xy-scatter
axes(handles.axScatter);
hS = scatter(xvals,yvals,5,cs);
xlabel(xtext,'FontSize',9);
ylabel(ytext,'FontSize',9);
hS.UserData = [xvals(:) yvals(:) roinums(:)];
hB = brush;
hB.Enable='on';
hB.ActionPostCallback = @brushingDone;
set(gcf,'UserData',handles);

function brushingDone(src,eventdata)
global beaddata beadGazerOpt beadSelection
ax = eventdata.Axes;
handles = get(gcf,'UserData'); % handles stored here
plt = ax.Children(1);
bdata = plt.BrushData;
lbox = handles.lbROIs;
% which axis are we on?
plttype = class(plt);
if contains(plttype,'Scatter')
    roinums = plt.UserData(:,3);
    beadSelection = roinums(logical(bdata));
    lbox.String = arrayfun(@(x) num2str(x), beadSelection,'UniformOutput',false);    
%     if handles.cbPosFilter.Value && ~isempty(beadGazerOpt.chosenPositions)
%         sel = beadGazerOpt.chosenPositions(bdata~=0);
%         lbox.String = arrayfun(@(x) num2str(x), sel,'UniformOutput',false);
%     else
%         lbox.String = arrayfun(@(x) num2str(x), find(bdata),'UniformOutput',false);
%     end
    if numel(lbox.String)>0
        lbox.Value = 1:numel(lbox.String);
        plotROIs(handles,false); % don't do data cursor mode
        handles.tRoiN.String = ['n = ' num2str(numel(lbox.String)) ...
            '/' num2str(numel(roinums)) ];
    else
        lbox.Value = [];
        cla(handles.axResponse);
    end
elseif contains(plttype,'Line')
    % find out which subplots have any points selected
    sel  = arrayfun(@(x) any(x.BrushData), ax.Children);
    sel  = flipud(sel); % because the plots are in reverse order
    % rois = arrayfun(@(x) x.UserData, ax.Children(sel));
    % disp(rois);
    dis  = lbox.Value; % all the plotted rois
    lbox.Value = dis(sel);  % the selected subset
    handles.tRoiN.String = ['nSel = '  num2str(sum(sel))];
    if ~isempty(lbox.Value)
        beadSelection = cellfun(@(s) str2double(s),lbox.String(lbox.Value));
        plotROIs(handles,false);
    else
        cla(handles.axResponse);
    end
end

function lbROIs_Callback(hObject, eventdata, handles)
handles.tRoiN.String = ['nSel = '  num2str(numel(hObject.Value))];
if ~isempty(hObject.Value), plotROIs(handles); end

function plotROIs(handles,enDataCursor)
global beaddata beadGazerOpt bg2

% when plotting values, valIdx says which column to use
valIdx = handles.lbColToPlot.Value;


str = handles.lbROIs.String;
sel = handles.lbROIs.Value;
if isempty(sel), return; end
rois = cellfun(@(x) str2double(x),str(sel));
try 
    sets = eval(handles.eSetsToPlot.String);
catch
    % disp('No sets specified - using all');
    sets = 1:size(beaddata.roiData,2);
end
cla(handles.axResponse);
data = beaddata.roiData;

% get the values/concentrations used for each set
manageConstraintsGUI(handles,valIdx,sets); 
sets = getSetConstraints(handles,sets);
nomValues =  beaddata.setInfo(sets,valIdx)';
handles.lbValues.String = arrayfun(@num2str,nomValues(:),'UniformOutput',false);

if handles.cbUseValues.Value  % cbPlotLogX
    xvals = nomValues;
    if handles.cbPlotLogX.Value
        xvals(xvals==0) = beadGazerOpt.logZero; % default value for "0" on log scale (value is raw, not log)
    end
else
    xvals = sets;
    % disp(nomValues(sets));
end

sw = isIntensity(handles);
if sw
    dvals = squeeze(data(rois,sets,2))';    % intensity
else
    dvals = squeeze(data(rois,sets,1))';    % lifetime
end

if handles.cbUseValues.Value 
    % plot the values in value order
    [xvals,idx] = sort(xvals);
    bg2.xvals = xvals;
    dvals = dvals(idx,:);
else
    bg2.xvals = nomValues;
end
if sw
    % calculate log10(F/F0)
    dvals = log10(dvals ./ repmat(dvals(1,:),size(dvals,1),1)); 
end  
bg2.yvals = dvals;

if handles.cbNoSymbols.Value
    p = plot(handles.axResponse,xvals, dvals); % plots each column as a line
else
    p = plot(handles.axResponse,xvals, dvals,'-o'); % plots each column as a line
end

% p=plot(xvals,squeeze(data(rois,sets,1))','-o'); % LT(,1) or intensity(,2)
for k=1:numel(p)
    p(k).UserData = rois(k);
end
if sw
    ylabel(handles.axResponse,'log10(F/F0)');
else
    ylabel(handles.axResponse,'Lifetime (ns)');
end

if isempty(sets), return; end

if handles.cbUseValues.Value 
    if handles.cbPlotLogX.Value
        % log axes
        set(handles.axResponse,'XScale','log');
        ax = axis(handles.axResponse);
        ax(1:2) = xvals([1 end]) .* [0.5 2];
        axis(handles.axResponse,ax);
    else
        ax = axis(handles.axResponse); % for later y-axis adjustment
        set(handles.axResponse,'XScale','linear');
    end
    set(handles.axResponse,'XTickMode','auto');
    set(handles.axResponse,'XMinorTick','on');
    xlabel(handles.axResponse,beaddata.valueNames{valIdx});
else
    % linear axes
    set(handles.axResponse,'XScale','linear');
    set(handles.axResponse,'XMinorTick','off')
    set(handles.axResponse,'XTickMode','manual');
    set(handles.axResponse,'XTick',sets);
    ax = axis(handles.axResponse);
    ax = ax + [-0.2 0.2 0 0];
    axis(handles.axResponse,ax);
    xlabel(handles.axResponse,'Set #');
end
if ax(4)-ax(3) < 0.4  % set a minimum extent for the y-axis
    ax(3) = ax(3)-0.1;
    ax(4) = ax(3)+0.5; 
    axis(handles.axResponse,ax);
end
if nargin>1, return; end  % avoids errors during brushing
dcm_obj = datacursormode(handles.figure1);
dcm_obj.UpdateFcn = {@myupdatefcn,rois};
% dcm_obj.Enable = 'on';
% datacursormode on;

function manageConstraintsGUI(handles,pltCol,sets)
global beaddata
bd = beaddata;
ncols = numel(bd.valueNames);
j = 1; % next constraint GUI to use
for k=1:ncols
    if k~=pltCol && k<=ncols
        js = ['Col' num2str(j)];
        hTxt = handles.(['t' js]);
        hPM  = handles.(['pm' js]);
        hTxt.String = bd.valueNames{k};
        hTxt.Visible = 'on';
        vals = bd.setInfo(sets,k);
        valstr = ['all values'; ...
            arrayfun(@(x) num2str(x) ,vals(:),'UniformOutput',false)];
        valstr = unique(valstr,'stable');
        hPM.String = valstr;
        if hPM.Value > numel(valstr), hPM.Value = 1; end
        hPM.Visible = 'on';
        hPM.UserData = k;
        j = j+1;
    end
end
% make the unused ones invisible
for k=j:4
    js = ['Col' num2str(k)];
    handles.(['t' js]).Visible = 'off';
    handles.(['pm' js]).Visible = 'off';
    handles.(['pm' js]).UserData = 0;
end
% set up the group by box
gblist = bd.valueNames;
gblist{pltCol} = ' - ';
handles.pmGroupBy.String = gblist;
handles.pmGroupBy.Value  = pltCol;

function sets = getSetConstraints(handles,setRange)
global beaddata bg2
bd = beaddata;
if nargin<2
    sets = true(size(bd.setInfo,1),1);
else
    sets = false(size(bd.setInfo,1),1);
    sets(setRange) = true;
end
% also save information about the condition to the global
cond = {'Plotting vs' bd.valueNames{handles.lbColToPlot.Value}};
for j=1:4
    js = ['Col' num2str(j)];
    hPM  = handles.(['pm' js]);
    col = hPM.UserData;
    if col>0 && hPM.Value>1   % if valid and not 'all'...
        selval = str2double(hPM.String{hPM.Value});
        sets = sets & bd.setInfo(:,col)==selval;
        cond{end+1,1} = bd.valueNames{col};
        cond{end,2}   = selval;
    end
end
sets = find(sets); % need a list of numbers
sets = sets(:)';
bg2.cond = cond;

function txt = myupdatefcn(~,event_obj,rois)
% Customizes text of data tips
pos = get(event_obj,'Position');
I = get(event_obj, 'DataIndex');
ud = event_obj.Target.UserData;
if numel(ud)==1
    roi = ud;
elseif isempty(ud) || size(ud,2)~=3
    roi = [];
else
    test = ud(:,1:2)==repmat(pos,size(ud,1),1);
    test = test(:,1) & test(:,2);
    roi = find(test);
end

txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['ROI: ',num2str(roi)]};


function cbLastMFirst_Callback(hObject, eventdata, handles)

% STILL TO DO:  
%   Multiple lines in response plotting e.g. sets = [1:2; 3:4]


function cbDTBasedOnVal_Callback(hObject, eventdata, handles)

function lbValues_Callback(hObject, eventdata, handles)

function lbColToPlot_Callback(hObject, eventdata, handles)
plotROIs(handles);

function cbPosFilter_Callback(hObject, eventdata, handles)
global beadGazerOpt
if hObject.Value && ~isfield(beadGazerOpt,'chosenPositions')
    pbChoosePos_Callback([],[],handles);
else
    pbReplot_Callback([],[],handles);
end

function pbChoosePos_Callback(hObject, eventdata, handles)
global beaddata beadGazerOpt
beadGazerOpt.chosenPositions = [];
if ~isempty(beaddata)
    beadGazerOpt.chosenPositions = brushXYselect;
end
pbReplot_Callback([],[],handles);

function cbNoSymbols_Callback(hObject, eventdata, handles)
plotROIs(handles);


function pmCol1_Callback(hObject, eventdata, handles)
plotROIs(handles);

function pmCol2_Callback(hObject, eventdata, handles)
plotROIs(handles);

function pmCol3_Callback(hObject, eventdata, handles)
plotROIs(handles);

function pmCol4_Callback(hObject, eventdata, handles)
plotROIs(handles);

function pbCopyData_Callback(hObject, eventdata, handles)
global bg2
% columns are X, Y1, Y2, ...
data = [bg2.xvals(:) bg2.yvals];
mat2clip(data);

function pbCopySummary_Callback(hObject, eventdata, handles)
global bg2
% columns are X, mean(Y), std(Y), N
data = [bg2.xvals(:) mean(bg2.yvals,2) std(bg2.yvals,0,2) repmat(size(bg2.yvals,2),numel(bg2.xvals),1)];
mat2clip(data);

function pbCopyCond_Callback(hObject, eventdata, handles)
global bg2
mat2clip(bg2.cond);


function pmGroupBy_Callback(hObject, eventdata, handles)
global beaddata bg2
bd = beaddata;
pltCol = handles.lbColToPlot.Value;
col2 = hObject.Value;
if pltCol==col2, return; end
% save the 'Use values' checkbox value, and set to true
useVals = handles.cbUseValues.Value;
handles.cbUseValues.Value = true;
% get the relevant values to scan (col2)
for conset=1:4
    hPM = handles.(['pmCol' num2str(conset)]);
    if hPM.UserData==col2
        break;
    end
end
legstr = hPM.String(2:end);  % for the legend
vals = arrayfun(@(x) str2double(x),legstr);
f = figure; hold on;
for k=1:numel(vals)
    hPM.Value = k+1; % make the selection
    % plot the current ROIs
    BeadGazer2('plotROIs',handles,false);
    errorbar(bg2.xvals, mean(bg2.yvals,2), ...
        std(bg2.yvals,0,2)/sqrt(size(bg2.yvals,2)));
end
leg = legend(legstr);
title(leg,hObject.String{col2});
xlabel(bd.valueNames{pltCol});
ylabel('Lifetime (ns)');
if handles.cbPlotLogX.Value
    set(gca,'XScale','log');
end
hObject.Value = pltCol;
% restore the 'Use values' checkbox value
handles.cbUseValues.Value = useVals;


function pmRoiGrpA_Callback(hObject, eventdata, handles)
performRoiGrpAction(hObject,eventdata,handles);
function pmRoiGrpB_Callback(hObject, eventdata, handles)
performRoiGrpAction(hObject,eventdata,handles);

function performRoiGrpAction(hObject, eventdata, handles)
global beadsurv
hList = handles.(['lb' hObject.Tag(3:end)]);
hInfo = handles.(['t'  hObject.Tag(3:end)]);
switch hObject.Value % which action to take
    case 2
        % Add plotted ROIs
        vals = hList.String;
        rois = handles.lbROIs.String;
        for k=handles.lbROIs.Value(:)'
            roi = str2double(rois(k));
            ij = bdTileIJforROI(roi);
            vals{end+1} = sprintf('%5i %3i  %3i', roi, ij(1), ij(2));
        end
        vals = unique(vals);
        hList.String = vals;
        % Also compile the list of unique tiles included
        tileI = cellfun(@(x) str2double(x(6:10)),vals);
        tileJ = cellfun(@(x) str2double(x(11:end)),vals);
        [tiles,~,ix] = unique([tileI(:) tileJ(:)],'rows','stable');
        hList.UserData.tiles = tiles;
        hList.UserData.freq  = accumarray(ix,1);
        hInfo.String = sprintf('%i ROIs in %i tiles',...
            numel(vals),size(tiles,1));
    case 3
        % Send to plotted
        vals = hList.String;
        %rois = cellfun(@(x) str2double(x(1:5),hList.String);
        rois = cellfun(@(x) x(1:5), vals, 'uni', false);
        handles.lbROIs.String = rois;
        handles.lbROIs.Value = 1:numel(rois);
    case 4
        % Make scope ctrl list
        % first we want to find the unique tiles
        tiles  = hList.UserData.tiles;
        oldfoldr = cd();
        try cd(beadsurv.filespec); end
        [fname,pname] = uigetfile('*.xml','experiment.xml file with full tile array','experiment.xml');
        if numel(fname)>1
            xmlstr = newXMLWithSubimages(fullfile(pname,fname),tiles);
            userroot = getenv('USERPROFILE');
            try
                cd('C:\Users\user\Documents');
            catch
                try cd([userroot '/Documents']); end
            end
            [fname,pname] = uiputfile('*.xml','New experiment.xml file for multiphoton  *** load afterwards in Capture Setup ***','experiment.xml');
            cd(oldfoldr);
            if numel(fname)>1  % nothing for cancel
                writecell(xmlstr,fullfile(pname,fname),'FileType','text');
                disp(['Wrote tile addresses for ' hObject.Tag(3:end) ' as ' fname ' in ' pname ]);
            end
        end
    case 5 
        % Clear the list...
        d = questdlg(['Ok to clear ' hObject.Tag(3:end) '?'],...
            'Confirm','Yes','No','No');
        switch d
            case 'Yes'
                hList.String = {};
                hInfo.String = '';
        end
    otherwise
        % nothing to do
end
hObject.Value = 1;

function data = getVals(sw,dataIn,sets,aggrType)
% dataIn = data(roi,set,datatype)
%   datatype defaults to lifetime unless sw is true
%   aggr type can be blank/absent or 'mean', 'min', 'max'
%   NB: intensities are usually ratioed over the first value
if sw, dtype=2; else dtype=1; end
data = dataIn(:,sets,dtype);
if nargin>=4
    switch aggrType
        case 'mean'
            data = mean(data,2);
        case 'max'
            data = max(data,[],2);
        case 'min'
            data = min(data,[],2);
    end
end

function tf = isIntensity(handles)
tf = handles.cbIntensity.Value;


function lbRoiGrpA_Callback(hObject, eventdata, handles)
function lbRoiGrpB_Callback(hObject, eventdata, handles)

function lbColToPlot_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function lbValues_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pmCol1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pmCol2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pmCol3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pmCol4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function lbROIs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eSetsToPlot_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eSetRange_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pmRoiGrpA_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pmRoiGrpB_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function lbRoiGrpB_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pmGroupBy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function lbRoiGrpA_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cbCriterion_Callback(hObject, eventdata, handles)
function cbAltX_Callback(hObject, eventdata, handles)
function cbAltY_Callback(hObject, eventdata, handles)
function eCriterion_Callback(hObject, eventdata, handles)
function eAltX_Callback(hObject, eventdata, handles)
function eAltY_Callback(hObject, eventdata, handles)


function showSetList(handles)
global beaddata
uit = handles.uitSetList;
nsets  = size(beaddata.setInfo,1);
ncols  = sum(cellfun(@(c) ~isempty(c), beaddata.valueNames));
cnames = beaddata.valueNames(1:ncols);
uit.ColumnName = cnames;
uit.Data = arrayfun(@(x) num2str(x), beaddata.setInfo(:,1:ncols), 'UniformOutput',false);
uit.RowName = arrayfun(@(x) ['s' num2str(x)], (1:nsets)', 'UniformOutput',false);


function pbPastiche_Callback(hObject, eventdata, handles)
global roiImages beadSelection
if isempty(roiImages)
    disp('No roiImages available from BeadSurveyorT (should be in global roiImages)');
    return
end
lbox = handles.lbROIs;
if ~isempty(lbox.Value)
    beadSelection = cellfun(@(s) str2double(s),lbox.String(lbox.Value));
    pas = beadPastiche(roiImages,beadSelection,4);
    figure(84); image(pas);
end


function pbRestrictImage_Callback(hObject, eventdata, handles)
% starting with the current selection of ROIs, choose only those with cs==1
global beaddata beadSelection
lbox = handles.lbROIs;
if ~isempty(lbox.Value) && isfield(beaddata,'cs')
    allROIs = cellfun(@(s) str2double(s),lbox.String);
    beadSelNow = allROIs(lbox.Value);
    beadSelection = beadSelNow(beaddata.cs(beadSelNow)==1);
    lbox.Value = lbox.Value(arrayfun(@(x) any(beadSelection==allROIs(x)),lbox.Value));
    handles.tRoiN.String = ['nSel = '  num2str(numel(lbox.Value))];
    if ~isempty(lbox.Value), plotROIs(handles); end
end


function pbInvertSelection_Callback(hObject, eventdata, handles)
global beaddata beadSelection
lbox = handles.lbROIs;
sel  = 1:numel(lbox.String);
sel(lbox.Value) = [];
lbox.Value = sel;
handles.tRoiN.String = ['nSel = '  num2str(numel(lbox.Value))];
if ~isempty(lbox.Value), plotROIs(handles); end


function cbIntensity_Callback(hObject, eventdata, handles)
