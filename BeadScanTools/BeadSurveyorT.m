function varargout = BeadSurveyorT(varargin)
% BEADSURVEYORT MATLAB code for BeadSurveyorT.fig  ThorimageLS-compatible version
%      BEADSURVEYORT, by itself, creates a new BEADSURVEYORT or raises the existing
%      singleton*.
%
%      H = BEADSURVEYORT returns the handle to a new BEADSURVEYORT or the handle to
%      the existing singleton*.
%
%      BEADSURVEYORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEADSURVEYORT.M with the given input arguments.
%
%      BEADSURVEYORT('Property','Value',...) creates a new BEADSURVEYORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BeadSurveyorT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BeadSurveyorT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BeadSurveyorT

% Last Modified by GUIDE v2.5 02-Mar-2020 10:04:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BeadSurveyorT_OpeningFcn, ...
                   'gui_OutputFcn',  @BeadSurveyorT_OutputFcn, ...
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


% --- Executes just before BeadSurveyorT is made visible.
function BeadSurveyorT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BeadSurveyorT (see VARARGIN)

% Choose default command line output for BeadSurveyorT
handles.output = hObject;

% addlistener(handles.uitTZero,'Data','PostSet',...
%     @(src,event) uitTZero_ListenerCallback(src,event,handles.uitTZero,handles));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BeadSurveyorT wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global xyzNow beadsurv bsGUI roiImages
evalin('base','global beadsurv bsGUI beaddata xyzNow dFLIMImages tExpt roiImages');
evalin('base','global beadImages_ref beadImages_now bsImages_now beadImages_pos');
xyzNow = [0 0 0];
handles.axes1.XTick=[];
handles.axes1.YTick=[];
handles.axes2.XTick=[];
handles.axes2.YTick=[];

bsGUI.handles = handles;  % save the GUI handles in a global
bsGUI.YDir2 = 'reverse';   % handles up-down reversal of secondary y-axis

if isempty(beadsurv)
    loadGlobalFromBSGUI(handles);
else 
    loadBSGUIFromGlobal(handles);
end
if ~isfield(beadsurv,'BL'), beadsurv.dZdXY = [0 0]; end

beadsurv.simulate = true;  % Thorimage: no direct moves used

function buttonDownCallback(src,~)
% TODO: Modify for communication with Thorimage
global bsGUI
tileIJ = src.UserData;
ok = questdlg(['Do you want to MOVE SCOPE to tile ' mat2str(tileIJ) ' ?']);
if strncmp(ok,'Yes',3)
    disp(['Moving to tile ' mat2str(tileIJ) ' ...']);
else
    return
end
xyzTarget = calcSubfieldPosition(tileIJ);
% relocate the outline box
[xs,ys] = patchFor(xyzTarget);
bsGUI.selector.Vertices = [xs(:) ys(:)];
uistack(bsGUI.selector,'top');
drawnow;
% use the CURRENT ACTUAL z-value
xyzNow = updateMotorPosition;
xyzTarget(3) = xyzNow(3);
launchMove(xyzTarget);
moveDone(xyzTarget);


function loadGlobalFromBSGUI(handles)
global beadsurv
% transfer values FROM the GUI
beadsurv.fov          = str2double(handles.eFOV.String);
beadsurv.ffCenters    = eval(handles.eFFCenters.String);
beadsurv.numSubfields = eval(handles.eNumSubfields.String);
beadsurv.offSubfield  = eval(handles.eSubfieldOffset.String);
beadsurv.fullField    = beadsurv.ffCenters + beadsurv.fov;
handles.eFullField.String = mat2str(beadsurv.fullField);
beadsurv.filespec       = handles.eFileSpec.String;
beadsurv.thresh       = str2double(handles.eThreshold1.String);
beadsurv.thresh2      = eval(handles.eThreshold2.String);
beadsurv.sens         = str2double(handles.eSens.String);
beadsurv.radrange     = eval(handles.eRadiusRange.String);
beadsurv.sourceChan   = eval(handles.eSourceCh.String);

function loadBSGUIFromGlobal(handles,partial)
global beadsurv beaddata
if nargin<2, partial=false; end
% transfer values TO the GUI
handles.eFOV.String            = num2str(beadsurv.fov);
handles.eFFCenters.String      = mat2str(beadsurv.ffCenters);
handles.eNumSubfields.String   = mat2str(beadsurv.numSubfields);
handles.eSubfieldOffset.String = mat2str(beadsurv.offSubfield);
beadsurv.fullField             = beadsurv.ffCenters + beadsurv.fov;
handles.eFullField.String      = mat2str(beadsurv.fullField);
if ~partial
    handles.eFileSpec.String       = beadsurv.filespec;
    handles.eThreshold1.String     = num2str(beadsurv.thresh);
    handles.eThreshold2.String     = mat2str(beadsurv.thresh2);
    handles.eSens.String           = num2str(beadsurv.sens);
    handles.eRadiusRange.String    = mat2str(beadsurv.radrange);
    handles.eSourceCh.String       = mat2str(beadsurv.sourceChan);
end
% updates the (grayed) corner positions
show_corners(handles);
if partial, return; end
% update the value names and current set values
try
    vnames = beaddata.valueNames;
    vdata = beaddata.setInfo;
    lastSet = size(vdata,1);
    handles.eSetNumber.String = num2str(lastSet+1);
    % construct and write the value tables
    for k=1:size(vnames,1)
        vnames{k,2} = vdata(lastSet,k);
    end
    handles.tblNamesValues.Data = vnames;
    handles.tblSetList.Data = vdata;
end

% --- Outputs from this function are returned to the command line.
function varargout = BeadSurveyorT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% Goto Corner pushbuttons
function pbTL_Callback(hObject, eventdata, handles)
global beadsurv xyzNow
launchMove(beadsurv.TL);
moveDone(beadsurv.TL);
% ask for focus at this position 
requestFocus;
beadsurv.TL(3) = xyzNow(3);
recalc_slopes;
show_corners(handles);

function pbTR_Callback(hObject, eventdata, handles)
global beadsurv xyzNow
launchMove(beadsurv.TR);
moveDone(beadsurv.TR);
% ask for focus at this position 
requestFocus;
beadsurv.TR(3) = xyzNow(3);
recalc_slopes;
show_corners(handles);

function pbBL_Callback(hObject, eventdata, handles)
global beadsurv xyzNow
launchMove(beadsurv.BL);
moveDone(beadsurv.BL);
% ask for focus at this position
requestFocus;
beadsurv.BL(3) = xyzNow(3);
recalc_slopes;
show_corners(handles);

function pbBR_Callback(hObject, eventdata, handles)
global beadsurv xyzNow
launchMove(beadsurv.BR);
moveDone(beadsurv.BR);
% ask for focus at this position
requestFocus;
beadsurv.BR(3) = xyzNow(3);
recalc_slopes;
show_corners(handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field size and position management  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eFFCenters_Callback(hObject, eventdata, handles)
global beadsurv
beadsurv.ffCenters = eval(handles.eFFCenters.String);
beadsurv.fullField = beadsurv.ffCenters + beadsurv.fov;

function eNumSubfields_Callback(hObject, eventdata, handles)
global beadsurv bsImages_now
% entered new number of subfields
cla(handles.axes1);
bsImages_now = {};
beadsurv.numSubfields = eval(handles.eNumSubfields.String);
if numel(beadsurv.numSubfields)~=2, disp('# subfields: specify as [12 24]'); return; end
beadsurv.ffCenters = (beadsurv.numSubfields-[1 1]) .* beadsurv.offSubfield;
loadBSGUIFromGlobal(handles);
recalc_corners(handles);
recalc_display(handles);

function recalc_slopes
% recalculates the z-offsets based on the focusing at the corners
global beadsurv
defined = isfield(beadsurv,{'BL' 'BR' 'TL' 'TR'});
if all(defined)
    dXYZ = mean([beadsurv.TL-beadsurv.BL; beadsurv.TR-beadsurv.BR]);
    dZdY = dXYZ(3)/dXYZ(2);
    dXYZ = mean([beadsurv.BR-beadsurv.BL; beadsurv.TR-beadsurv.TL]);
    dZdX = dXYZ(3)/dXYZ(1);
else
    dZdY = 0;
    dZdX = 0;
end
beadsurv.dZdXY = [dZdX dZdY];

function recalc_corners(handles)
global beadsurv
% recalculate the corners relative to BL, 
%     based on subfield offset and number of subfields
% interpolate Z values based on old values
extent  = (beadsurv.numSubfields - [1 1]) .* beadsurv.offSubfield;
corners = {'BL' 'BR' 'TL' 'TR'};
defined = isfield(beadsurv,corners);
recalc_slopes;
dZdX = beadsurv.dZdXY(1);
dZdY = beadsurv.dZdXY(2);
if defined(1)  % we do have a bottom left to use
    beadsurv.BR = beadsurv.BL + [extent(1) 0 dZdX*extent(1)];
    beadsurv.TL = beadsurv.BL + [0 -extent(2) -dZdY*extent(2)];
    beadsurv.TR = beadsurv.BL + [extent.*[1 -1] extent*[dZdX; -dZdY]];
end
show_corners(handles);

function show_corners(handles)
global beadsurv
corners = {'BL' 'BR' 'TL' 'TR'};
defined = isfield(beadsurv,corners);
for k=1:4
    if defined(k)
        handles.(['eXYZ_' corners{k}]).String = ...
            mat2str(round(beadsurv.(corners{k}),1));
    end
end

function eSubfieldOffset_Callback(hObject, eventdata, handles)
global beadsurv
subfOffset = eval(handles.eSubfieldOffset.String);
if numel(subfOffset)~=2, disp('subfield offset: specify as [offX offY]'); return; end
beadsurv.offSubfield = subfOffset; 
fullWH = (beadsurv.numSubfields - [1 1]) .* subfOffset;
beadsurv.ffCenters = fullWH;
beadsurv.fullField = fullWH + beadsurv.fov;
loadBSGUIFromGlobal(handles);
recalc_corners(handles);
recalc_display(handles)

function eFOV_Callback(hObject, eventdata, handles)
global beadsurv
beadsurv.fov = str2double(hObject.String);
% recalc the full field, but not the position of the corners
loadBSGUIFromGlobal(handles);
recalc_display(handles)



function pbTestSegment_Callback(hObject, eventdata, handles)
global tExpt bsGUI
% start by making sure we have an open tExpt
fpath    = handles.eFileSpec.String;  % SHOULD BE A PATH NOW
if isempty(tExpt) || isempty(tExpt.info) || ~strcmpi(tExpt.exptPath,fpath)
    if fpath(end)~='\', fpath = [fpath '\']; end
    tExpt    = TIExpt(fpath);
end
tInfo    = tExpt.info;
if isempty(tInfo), dispr('Data unavailable!'); beep(); return; end

nch = numel(tExpt.info.channels);
chan = eval(handles.eSourceCh.String);
[~,idx]   = tExpt.serialFromTileIJ([str2double(handles.eTestCol.String),...
    str2double(handles.eTestRow.String)]);
[dfh,dfi] = tExpt.getData(idx);
processTZero(dfh,dfi,chan,handles);
LUTArray = getLUTValues(handles); % [1.8 2.7; 0 5];
greyLUT = LUTArray(2,:);
if handles.cbLT.Value
    if ~handles.cbSmooth.Value
        cImg = dfi(chan(1)).LifetimeRGBImage(LUTArray,0);
    else
        sm1 = str2double(handles.eSmooth.String);
        sm2 = max(3,min(9,1+2*round(sm1*1.3)));
        %disp(['Smooth ' num2str(sm1) ' ' num2str(sm2)]);
        cImg = dfi(chan(1)).LifetimeRGBImage(LUTArray,0,sm1,sm2);
    end
else
    cImg = dfi(chan(1)).CDataS(greyLUT(1),greyLUT(2));
end

%622 greyLUT = [0 30];

% analyse for ROIs, report on the number
[centers,radii,inten,ltime,nphot,cOs] = ...
    segmentThisImage(dfi(chan(1)),cImg,handles); % TODO: param1=intensity image, param2=RGB-lifetime-image
nROIs = numel(radii);
disp(['nROIs = ' num2str(nROIs)]);



% segments the current image
function pbRedoSubfield_Callback(hObject, eventdata, handles)
global dfiNow dFLIMImageDisplay beadsurv
chan     = eval(handles.eSourceCh.String);
beadsurv.sourceChan = chan;
segmentThisImage(dfiNow(chan(1)),dFLIMImageDisplay{chan(1),3},handles);

function [ctrs,radii,inten,ltime,nsingl,roiImage] = segmentThisImage(dfiNow,dispImage,handles)
% given the current subfield image, recalculate the circular ROIs and display
global beadsurv bsGUI
% changed to NMG (gy 20191210)
%imgNow   = flipud(dfiNow.nmg); % FLIPPED because of axis orientation chg; was dfiNow.img;
imgNow   = (dfiNow.nmg); % NOT FLIPPED 
nFrames  = dfiNow.nFrames;
radrange = eval(handles.eRadiusRange.String);
if handles.cbRadiusShrink.Value
    beadsurv.radshrink = str2double(handles.eRadiusShrink.String);
else
    beadsurv.radshrink = 0;
end
sens     = str2double(handles.eSens.String);
% first threshold number is minimum for drawing
thresh   = str2double(handles.eThreshold1.String);
% second threshold field specifies which to select for numbering/data coll
%  either a single number (min value) or two numbers [min max]
thresh2  = eval(handles.eThreshold2.String); 
beadsurv.thresh  = thresh;
beadsurv.thresh2 = thresh2;
beadsurv.sens    = sens;
beadsurv.radrange= radrange;
%tZero = 2.84;
s=warning('off','images:imfindcircles:warnForSmallRadius');
[ctrs,radii] = imfindcircles(imgNow>thresh,radrange,'Sensitivity',sens,'method','twostage');
warning(s);

% apply the ROI shrink to avoid the shell fluorescence (if checked)
radii = radii - beadsurv.radshrink;

axes(handles.axes2); cla; 
imshow((dispImage));  % not flipud
set(handles.axes2,'YDir',bsGUI.YDir2); % for Thorimage
handles.axes2.XTick=[];
handles.axes2.YTick=[];
%imshow(uint8(255*mat2gray(imgNow))); % imshow(imgNow>thresh); 

[xx,yy] = meshgrid(1:size(imgNow,2),1:size(imgNow,1));

if numel(radii)==0, inten=[]; ltime=[]; nsingl=[]; roiImage={}; return; end

lt=zeros(size(radii));
iten = lt;
nph = lt;
% core/shell, 1/2/3/4 moments of intensity distribution among the pixels
roiImage = cell(size(radii));

for k=1:numel(radii)
    % can improve speed here using calcBeadSurveyorMask; just need meshgrid
    % xxL and yyL precalculated
    maskFull = hypot(xx-ctrs(k,1),yy-ctrs(k,2))<radii(k);
    [inten,ltime,nphotFull,nsingl] = dfiNow.roiValues(maskFull);
    lt(k) = ltime;
    iten(k) = inten;
    nph(k) = nsingl;
%     maskCore = hypot(xx-ctrs(k,1),yy-ctrs(k,2))<0.5*radii(k); % middle 25% of area
%     [~,~,nphotCore] = dfiNow.roiValues(maskCore);
%     % calculate coreOverShell mean intensities
%     cOs(k,1) = (nphotCore/sum(maskCore(:))) / ...
%         ((nphotFull - nphotCore)/(sum(maskFull(:))-sum(maskCore(:))));
    % calculate non-uniformity (sd/mean)
%     intensities = dfiNow.img(maskFull);
%     cOs(k,2) = inten;
%     cOs(k,3:5) = [std(intensities(:)) skewness(intensities(:)) kurtosis(intensities(:))] ...
%         / inten;
    %
    % prepare an image excerpt for each ROI
    %
    approxCtr = round(ctrs(k,:));
    spread    = ceil(radii(k)+1);
    xy = [max(1,approxCtr - [spread spread]) ...
        min([size(imgNow,2) size(imgNow,1)], approxCtr + [spread spread])];
    roiImg    = dfiNow.img(xy(2):xy(4),xy(1):xy(3));
    roiCtr    = (min(1,approxCtr-spread*[1 1])+spread)+ctrs(k,:)-approxCtr;
    roiImage{k} = {roiImg roiCtr radii(k)};
end
% toc,
if numel(thresh2)==2
    %sel = (nph > thresh2(1)/nFrames) & (nph <= thresh2(2)/nFrames);
    sel = (nph > thresh2(1)) & (nph <= thresh2(2));
else
    %sel = nph > thresh2/nFrames;
    sel = nph > thresh2;
end
ctrs  = ctrs(sel,:);
radii = radii(sel);
inten = iten(sel);
ltime = lt(sel);
nsingl = nph(sel);
roiImage = roiImage(sel);

if handles.cbLT.Value
    viscircles(ctrs,radii,'EdgeColor','w','LineWidth',0.5);
else
    viscircles(ctrs,radii,'EdgeColor','r','LineWidth',0.5);
end
if handles.cbShowID.Value
    for k=1:numel(radii)
        %if iten(k) > thresh2/nFrames
        text(1+ctrs(k,1)-0.75*radii(k),ctrs(k,2), sprintf('%2i',k),...
            'Color','g','UserData',[k inten(k) ltime(k) nsingl(k)],...
            'ButtonDownFcn',@roiButtonDown);
        % viscircles(ctrs(k),radii(k));
        %end
    end
end

% lt = lt(iten>thresh*2/nFrames);
% nph = nph(iten>thresh*2/nFrames);
% iten = iten(iten>thresh*2/nFrames);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image set collection commands       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roiButtonDown(src,ev)
ud   = src.UserData;
hAx  = src.Parent;
pos  = src.Position(1:2)';
pos(2) = hAx.YLim(2) - pos(2);
% lims = [hAx.XLim; hAx.YLim];
% posR = (pos-lims(:,1))./(lims(:,2)-lims(:,1));
% posR(2) = 1-posR(2);
hh = annotation('textbox',axxy2figxy(hAx,[pos' .2 .1]),'FitBoxToText','on',...
    'BackgroundColor','w','FontSize',8,...
    'String',{num2str(ud(1)); ['Intensity ',num2str(ud(2),3)];...
       ['Lifetime ' num2str(ud(3),3)]; ['nPhotons ' num2str(ud(4))]});
pause(2);
delete(hh);

function pbCollectSet_Callback(hObject, eventdata, handles)
global beadsurv tExpt bsGUI
% start by monitoring files
fpath    = handles.eFileSpec.String;  % SHOULD BE A PATH NOW
if fpath(end)~='\', fpath = [fpath '\']; end
tExpt    = TIExpt(fpath);
tInfo    = tExpt.info;
if isempty(tInfo), dispr('Data unavailable!'); beep(); return; end
if ~isfield(tInfo,'tiles'), dispr('Not a tiled dataset'); beep(); return; end
%%%%%%%%%%%% monitorImageFiles(fspec);
% now get and read the experiment config from Thorimage
% i = strfind(fspec,'\');
% fpath = fspec(1:i(end));
% tExpt = readThorimageExperimentFile([fpath 'Experiment.xml']);
bsGUI.pathspec = fpath;  %TODO: CHECK

setNum = str2double(handles.eSetNumber.String);
if setNum==1
    % interpret the tile layout
    beadsurv.fov          = tInfo.fieldsize(1);
    beadsurv.numSubfields = tInfo.tilesXY;
    beadsurv.offSubfield  = abs(tInfo.tileDxy);
    beadsurv.ffCenters    = (tInfo.tilesXY-1) .* [1 -1] .* tInfo.tileDxy;
    beadsurv.TL           = tInfo.tile0xy;
    beadsurv.TR           = beadsurv.TL + [1 0] .* beadsurv.ffCenters;
    beadsurv.BL           = beadsurv.TL + [0 -1] .* beadsurv.ffCenters;
    beadsurv.BR           = beadsurv.TL + [1 -1] .* beadsurv.ffCenters;
end
% display the info
loadBSGUIFromGlobal(handles,true);
prepForSet(setNum);
% monitor the files and deliver them for processing
tExpt.deliverTileData(@processBeadImage,setNum);

% importBeadImages(fspec,setNum,@processBeadImage);

% autoincrement set number and target path
handles.eSetNumber.String = num2str(setNum+1);
if fpath(end-4)=='_' % there is already a number at the end of the path 
    oldPathNumber = str2double(fpath(end-3:end-1));
    newPathNumber = oldPathNumber+1;
    fpath(end-3:end-1) = num2str(newPathNumber,'%03u');
    handles.eFileSpec.String = fpath;
else
    fpath(end:end+4) = '_001\';
    handles.eFileSpec.String = fpath;    
end
    
function prepForSet(setNum)
global beadsurv bsGUI beaddata
global beadImages_ref beadImages_now bsImages_now beadImages_pos
handles = bsGUI.handles;
setInfo = handles.tblNamesValues.Data;
if setNum==1
    % starting fresh:  reset the image records
    beadImages_ref = {};
    beadImages_now = {};
    beadImages_pos = {};
    bsImages_now = {};
    beaddata = [];
    beaddata.setInfo = cell2mat(setInfo(:,2))';
    handles.tblSetList.Data = beaddata.setInfo;
    beaddata.valueNames = setInfo(:,1);
    handles.tblSetList.ColumnName = beaddata.valueNames;
else
    beaddata.setInfo(setNum,:) = cell2mat(setInfo(:,2))';
    handles.tblSetList.Data = beaddata.setInfo;
end
% set up the display axes
%   units are um, need to figure out XDir and YDir for scope movement
axes(handles.axes1);
recalc_axlimits(handles);
xyTarget = calcSubfieldPosition([1 1]);
% create the blue outline box for the currently acquired image, if needed
if ~isfield(bsGUI,'selector') || ~isvalid(bsGUI.selector)
    % prevent imshow from deleting everything!
    iptsetpref('ImshowAxesVisible','on');
    [xs,ys] = patchFor(xyTarget);
    bsGUI.selector = patch(xs,ys,'red','FaceColor','none','EdgeColor','red');
    hold on;
end

function processBeadImage(tieObj, idx, setNum)
% called from importBeadImages, with the tile image filenames in sequence
% note that idx(2) is the serial number of the tile
global beadsurv bsGUI beaddata tExpt roiImages
global beadImages_ref beadImages_now bsImages_now beadImages_pos
handles = bsGUI.handles;

% retrieve the data from the newly available file
[dfh,dfi] = tieObj.getData(idx);
ijNow     = tieObj.tileIJ(idx);

% calculate the xy display position for this tile
xyThis = calcSubfieldPosition(ijNow);

% calculate the index of this tile, in order of acquisition
%  NB - this is based on the original full tile array, even when
%       subsequent sets use only limited single tiles
iTile = tileAcqorderFromIJ(ijNow,beadsurv.numSubfields);

fov2   = beadsurv.fov/2;
chan = eval(handles.eSourceCh.String);
beadsurv.sourceChan = chan;



% function processBeadImage(imfname,ijNow,info,setNum)
% % called from importBeadImages, with the tile image filenames in sequence
% % note that info(3) is the serial number of the tile
% global beadsurv bsGUI beaddata tExpt
% global beadImages_ref beadImages_now bsImages_now beadImages_pos
% handles = bsGUI.handles;
% disp(['Processing file: ' imfname]);
% cletters = chanLetters(chan);

imgdim  = tieObj.info.pixelsXY([2 1]); % [dFLIMImages(chan(1)).nLines dFLIMImages(chan(1)).nPixels];
imgOrig = [-1 -1]*fov2;  % position of lowest image pixel, in um
pxScale = 2*fov2 ./ imgdim;  % pixel size in um

% relocate the outline box
[xs,ys] = patchFor(xyThis);
bsGUI.selector.Vertices = [xs(:) ys(:)];
uistack(bsGUI.selector,'top');

%% image retrieval

% % get the file images:
% % [dfi,dfh] = loadRaw_Thorimage(imfname,imgsz,channums,nAveraged,tPeakIRFs);
% [dfi,dfh] = loadRaw_Thorimage(imfname,tExpt.pixelsXY,tExpt.channels,tExpt.averageNum,repmat(0.8,1,numel(tExpt.channels)));

nch = numel(tieObj.info.channels);

% % TODO:  This is the non-FLIM version
% images = cell(numel(chan));
% for iCh = 1:numel(chan)
%     fname    = imfname;
%     fname(5) = cletters(iCh);
%     fullname = [bsGUI.pathspec fname];
%     iFile    = Tiff(fullname,'r');
%     images{iCh} = iFile.read;  % reads an image from the Tiff file
% end

beadImages_now{ijNow(1),ijNow(2)} = dfi;  % save the current image(s)
beadImages_pos{ijNow(1),ijNow(2)} = [xyThis iTile];

% erase the old image at this position
%  (more efficient alternative is to replace its data)
try
    delete(bsImages_old{ijNow(1),ijNow(2)});
catch ME
    % disp(ME);
end
    
%% image display
% % TODO: for now just using a bw image
% cImg = double(images{1});
% cImg = cImg - min(cImg(:));
% cImg = uint8(256*cImg/max(cImg(:)));
% cImg = flipud(ind2rgb(cImg,bsGUI.handles.figure1.Colormap));

LUTArray = getLUTValues(handles);
processTZero(dfh,dfi,chan,handles);

%cImg = flipud(dfi(chan(1)).LifetimeRGBImage(LUTArray));
if ~handles.cbSmooth.Value
    cImg = dfi(chan(1)).LifetimeRGBImage(LUTArray,0);
else
    sm1 = str2double(handles.eSmooth.String);
    sm2 = max(3,min(9,1+2*round(sm1*1.3)));
    %disp(['Smooth ' num2str(sm1) ' ' num2str(sm2)]);
    cImg = dfi(chan(1)).LifetimeRGBImage(LUTArray,0,sm1,sm2);
end
%cImg = flipud(cImg);

if setNum==1
    % save the reference image
    beadImages_ref{ijNow(1),ijNow(2)} = dfi;
    
    % note that this imshow reverses the YDir of the display!
    %   but it may not when we use the 'CData' spec
    h = image(handles.axes1,'CData',flipud(cImg),'XData',xyThis(1)+[-fov2 fov2], ...
        'YData',xyThis(2)+[-fov2 fov2]);
    h.AlphaData = 1 - 0.5*handles.cbVisOverlap.Value;
    h.UserData = ijNow;
    h.ButtonDownFcn = @buttonDownCallback;
    try delete(bsImages_now{ijNow(1),ijNow(2)}); end
    bsImages_now{ijNow(1),ijNow(2)} = h;
    % save for quick reference:  
    beadImages_pos{ijNow(1),ijNow(2)} = [xyThis iTile];
else
    % just replace the data of the image object
    bsImages_now{ijNow(1),ijNow(2)}.CData = flipud(cImg);
end
recalc_axlimits(handles);
drawnow;

%% image processing
%622 greyLUT = [0 30];
LUTArray = getLUTValues(handles); % [1.8 2.7; 0 5];
greyLUT = LUTArray(2,:);

if setNum==1  % THE FIRST SET BEING COLLECTED - NEED TO SEGMENT
    % analyse for ROIs, report on the number
    % adding roiImg{k} = {imgArray,ctr,radius}
    if handles.cbLT.Value
        [centers,radii,inten,ltime,nphot,roiImg] = ...
            segmentThisImage(dfi(chan(1)),dfi(chan(1)).CDataS(greyLUT(1),greyLUT(2)),handles); % TODO: param1=intensity image, param2=RGB-lifetime-image
    else
        [centers,radii,inten,ltime,nphot,roiImg] = ...
            segmentThisImage(dfi(chan(1)),cImg,handles); % TODO: param1=intensity image, param2=RGB-lifetime-image
    end
    nROIs = numel(radii);
    if iTile==1
        firstROI = 1;
        roiImages = {};
    else
        firstROI = sum(beaddata.tileInfo(iTile-1,1:2));
    end
    %
    % .tileInfo has info on each tile and its ROIs
    %
    beaddata.tileInfo(iTile,1) = firstROI;
    beaddata.tileInfo(iTile,2) = nROIs;
    beaddata.tileInfo(iTile,3:4) = ijNow;
    beaddata.tileInfo(iTile,5:7) = [xyThis 0];  %TODO CHECK using 0 for z?
    handles.uitable3.Data = {iTile,ijNow(1),ijNow(2),firstROI,nROIs,firstROI+nROIs-1};
    rois = firstROI:(firstROI+nROIs-1);
    if nROIs>0
        %
        % store both the relative and absolute positions
        %   in .roiInfo(roi#,:) = [x y radius xAbs yAbs];
        %
        beaddata.roiInfo(rois, 1:5) = [centers radii ...
            (repmat(pxScale,nROIs,1) .* centers) + ...
            repmat(xyThis+imgOrig,nROIs,1)  ];
        %
        % .roiData(roi#,set#,:) = [lifetime, intensity, nPhotons]
        %
        beaddata.roiData(rois,setNum,1:3)=[ltime(:) inten(:) nphot(:)];
        if numel(chan)>1 % get data points for additional channels
            [xxL,yyL] = meshgrid(-10:10,-10:10);  % local grid; center +/- 10 pixels
            szimg = fliplr(size(dfi(chan(1)).img));
            for iROI = rois
                xyr = beaddata.roiInfo(iROI,:);
                mask0 = calcBeadSurveyorMask(xyr,[0 0],xxL,yyL,szimg);
                for iChan = 2:numel(chan)
                    % ********* NEEDS WORK **********
                    [inten,ltime,~,nphot] = dfi(chan(iChan)).roiValues(mask0);
                    beaddata.roiDataExt(iROI,setNum,1:3,iChan-1) = [ltime inten nphot];
                end
            end
        end
        %
        % .coreOverShell(roi#) = coreOverShell for primary channel
        %
        roiImages(rois) = roiImg;
    end
else
    % NON-INITIAL SET:  register the image and analyze ROIs
    % adjust the image to the reference image
    totalShift = -findPixelShift(...
        beadImages_ref{ijNow(1),ijNow(2)}(chan(1)).img, ...
        beadImages_now{ijNow(1),ijNow(2)}(chan(1)).img,5);
    disp(['>>>>>>>>shft ' mat2str(totalShift)]);
    dxy = totalShift; % [0 0];  % TODO - might need to be flipped
    % analyze the previously discovered ROIs
    firstROI = beaddata.tileInfo(iTile,1);
    nROIs = beaddata.tileInfo(iTile,2);
    roiRange = firstROI:(firstROI+nROIs-1);
    handles.uitable3.Data(1:5) = {iTile,ijNow(1),ijNow(2),firstROI,nROIs};
    
    % display the new image with the ROIs circled
    axes(handles.axes2); cla;
    imshow(dfi(chan(1)).CDataS(greyLUT(1),greyLUT(2)));
    %imshow(uint8(255*mat2gray(dFLIMImages(chan).img)));
    viscircles(beaddata.roiInfo(roiRange,1:2)+repmat(dxy,numel(roiRange),1),beaddata.roiInfo(roiRange,3),...
            'LineWidth',0.5,'EdgeColor','b');
    set(handles.axes2,'YDir',bsGUI.YDir2); % for Thorimage
    handles.axes2.YTick=[];
    
    % calculate ROI masks using only a subarray for the circle calculation
    [xxL,yyL] = meshgrid(-10:10,-10:10);  % local grid; center +/- 10 pixels
    szimg = fliplr(size(dfi(chan(1)).img));
    for iROI = roiRange
        xyr = beaddata.roiInfo(iROI,:);
        mask0 = calcBeadSurveyorMask(xyr,dxy,xxL,yyL,szimg);
        [inten,ltime,~,nphot] = dfi(chan(1)).roiValues(mask0);
        beaddata.roiData(iROI,setNum,1:3)=[ltime inten nphot];
        % define additional data structure for additional channels
        for iChan = 2:numel(chan)
            [inten,ltime,~,nphot] = dfi(chan(iChan)).roiValues(mask0);
            beaddata.roiDataExt(iROI,setNum,1:3,iChan-1) = [ltime inten nphot];
        end
    end
end

% save if we're done
% last tile case is different if we have a full array or not
if handles.cbAutosaveBeaddata.Value==1 && ( ...
        (isfield(tExpt.info,'tilesXY') && iTile==prod(beadsurv.numSubfields)) ...
        || (~isfield(tExpt.info,'tilesXY') && idx(2)==numel(tExpt.info.tiles)))
    % we save beadsurv and beaddata globals
    % new file at the end of each set (even though redundant)
    fname  = num2str(setNum,'_beaddata_s%02i');
    fname2 = 'roiImages';
    try
        fname  = fullfile(bsGUI.pathspec, fname);
        fname2 = fullfile(bsGUI.pathspec, fname2);
    catch
        % just save in the current folder if no savepath
    end
    % for first set, save the reference roiImages and the full reference
    % images
    if setNum==1
        save(fname,'beaddata','beadsurv','beadImages_ref');
        save(fname2,'roiImages','roiImages');
    else
        % in the future, omit the reference images for later sets?
        save(fname,'beaddata','beadsurv','beadImages_ref');
    end
end

    
function cletters = chanLetters(chan)
cl = 'ABCD';
cletters = cl(chan);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Master scanning function - gets image, segments (or registers  %
%   to the ref image), and analyzes ROI values                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 201705: right now only the 'Source channel' is analyzed
%

function [xs ys] = patchFor(xyzTarget)
global beadsurv
off2 = beadsurv.fov/2;
pos = xyzTarget(1:2);% - beadsurv.TL(1:2);
% pos is the position relative to the top left
% this is the center of the image, but the lower left corner of the patch
corners = pos + [1 1] * off2 ;  % first row
corners(2,:) = pos + [1 -1] * off2;
corners(3,:) = pos + [-1 -1] * off2;
corners(4,:) = pos + [-1 1] * off2;
xs = corners(:,1)';
ys = corners(:,2)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Tile position utilities        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xyzTarget = moveToTile(tileIJ)
xyzTarget = calcSubfieldPosition(tileIJ);
launchMove(xyzTarget);

function ijStart = startingIJ
global beadsurv
ijStart = [1 beadsurv.numSubfields(2)];

function xyz = calcSubfieldPosition(subIJ)
% subfield IJs are measured relative to Top Left
global beadsurv
bs = beadsurv;
% sign reversal for Thorimage setup
offset = [1 -1] .* beadsurv.offSubfield .* (subIJ - [1 1]);
%TODO zOffset = offset * beadsurv.dZdXY(:);
xyz = bs.TL + offset; %TODO + [offset zOffset];


function subIdx = idxSubfield(subIJ)
% return a linear index from BL=1, BR=nSubX, adding nSubX for each row.
%   (matrix indexed from 1 at bottom left)
global beadsurv
nSubf = beadsurv.numSubfields;

% get the row offset first
rowrel    = nSubf(2) - subIJ(2);
rowrising = 2*floor(rowrel/2)==rowrel;

if rowrising
    subIdx = rowrel*nSubf(1) + subIJ(1);
else
    subIdx = rowrel*nSubf(1) + 1 + nSubf(1) - subIJ(1);
end

function cbVisOverlap_Callback(hObject, eventdata, handles)
% visualize the image overlaps by making them transparent
global bsImages_now
if hObject.Value, alpha=0.5; else alpha=1; end
for k=1:numel(bsImages_now)
    h = bsImages_now{k};
    if ishandle(h), h.AlphaData=alpha; end
end

function recalc_display(handles)
% redisplays all the current images
global bsImages_now beadsurv
fov2 = beadsurv.fov/2;
axes(handles.axes1);
recalc_axlimits(handles);
for i=1:size(bsImages_now,1)
    for j=1:size(bsImages_now,2)
        h = bsImages_now{i,j};
        if ishandle(h)
            xyCenter = calcSubfieldPosition([i j]);
            h.XData = xyCenter(1) + [-fov2 fov2];
            h.YData = xyCenter(2) + [-fov2 fov2];
        end
    end
end

function recalc_axlimits(handles)
% changed for Thorimage version
%  because stage directions run positive toward top left
% use a NON-reversed axis
% then we will need to flip images ud to display correctly
global beadsurv bsGUI
axes(handles.axes1);
axis manual;
% set(handles.axes1,'YDir','reverse'); % works best for image displays
fov2 = beadsurv.fov/2;
bleft = beadsurv.BL;
axmax = beadsurv.ffCenters + fov2;
axlim = [-fov2 axmax(1) -fov2 axmax(2)] + ...
    [bleft(1) bleft(1) bleft(2) bleft(2)];
axis(axlim);    
set(handles.axes1,'YDir','normal'); % for Thorimage
set(handles.axes2,'YDir',bsGUI.YDir2); % for Thorimage
handles.axes2.YTick=[];

function cbAbortCollection_Callback(hObject, eventdata, handles)
global tExpt
if hObject.Value
    disp('Abort requested after timeout...');
    tExpt.abort = true;
    pause(5);
    hObject.Value = false;
end


% MICROSCOPE ACTION FUNCTIONS

function launchMove(xyz)
global beadsurv
if beadsurv.simulate
    disp(['Launch move to ' mat2str(xyz)]);
else
    xyzStartMove(xyz);
end

function ok = moveDone(xyz)
global xyzNow beadsurv state
if beadsurv.simulate
    pause(0.3);
    disp(['Finished move to ' mat2str(xyz)]);
    xyzNow = xyz;
else
    xyzFinishMove(true,[3 3 3],30);  % timeout of 30 is about 9 seconds
    xyzNow = [state.motor.absXPosition state.motor.absYPosition state.motor.absZPosition];
end
ok = 1;

function updateXYZ
global beadsurv state xyzNow
if ~beadsurv.simulate
    updateMotorPosition;
    xyzNow = [state.motor.absXPosition state.motor.absYPosition state.motor.absZPosition];
end

function eFileSpec_Callback(hObject, eventdata, handles)
global beadsurv
beadsurv.filespec = hObject.String;

function pbSaveBeadData_Callback(hObject, eventdata, handles)
global beadsurv beaddata beadImages_ref
try 
    fname=fullfile(state.files.savePath, [beadsurv.filespec '_beaddata']);
catch
    fname=[beadsurv.filespec '_beaddata'];
end
[fname pname]=uiputfile('*.mat','Save bead data',fname);
if ~isequal(fname,0) && ~isequal(pname,0)
    save(fullfile(pname,fname),'beaddata','beadsurv','beadImages_ref');
end

function pbShowROI_Callback(hObject, eventdata, handles)
global beadsurv beaddata beadImages_ref beadImages_now bsGUI
LUTArray = getLUTValues(handles); % [1.8 2.7; 0 5];
greyLUT = LUTArray(2,:);
chan = eval(handles.eSourceCh.String);
answer = inputdlg('ROI number to show?','BeadSurveyor:ShowROI',1);
if ~isempty(answer)
    iROI = str2double(answer);
    iTile = find(beaddata.tileInfo(:,1)<=iROI,1,'last');
    if isempty(iTile), error('BeadSurveyor:ShowROI - ROI number out of range'); return; end
    ij = beaddata.tileInfo(iTile,3:4);
    imgOrig = [-1 -1]*beadsurv.fov/2;  % position of lowest image pixel, in um
    roiInfo = beaddata.roiInfo(iROI, 1:5);
    disp(['ROI# ' num2str(iROI) ' idx ' num2str(iTile) ' ij ' mat2str(ij)]);
    disp(mat2str(roiInfo));
    % = [centers radii ...
%                 (repmat(pxScale,nROIs,1) .* centers) + ...
%                 repmat(xyzThis(1:2)+imgOrig,nROIs,1)  ];
    dfi =  beadImages_now{ij(1),ij(2)}; %622 img = beadImages_now{ij(1),ij(2)}.img;
    axes(handles.axes2); cla; imshow(dfi(chan(1)).CDataS(greyLUT(1),greyLUT(2))); %622 imshow(uint8(255*mat2gray(img)));
    set(handles.axes2,'YDir',bsGUI.YDir2); % for Thorimage
    handles.axes2.YTick=[];
    viscircles(roiInfo(1:2),roiInfo(3),'LineWidth',1,'EdgeColor','m');
    disp(['xy-center(frame) =' mat2str(round(roiInfo(4:5),1))]);
    disp(['xy-center(roi)   =' mat2str(round(roiInfo(1:2)+roiInfo(4:5),1))]);
end

function pbGoToROI_Callback(hObject, eventdata, handles)
% TODO: Modify to deliver XY to Thorimage and request move
global beadsurv beaddata beadImages_ref beadImages_now bsGUI
LUTArray = getLUTValues(handles); % [1.8 2.7; 0 5];
greyLUT = LUTArray(2,:);
chan = eval(handles.eSourceCh.String);
answer = inputdlg('ROI number to go-to?','BeadSurveyor:GoToROI',1);
if ~isempty(answer)
    iROI = str2double(answer);
    iTile = find(beaddata.tileInfo(:,1)<=iROI,1,'last');
    if isempty(iTile), error('BeadSurveyor:GoToROI - ROI number out of range'); return; end
    ij = beaddata.tileInfo(iTile,3:4);
    imgOrig = [-1 -1]*beadsurv.fov/2;  % position of lowest image pixel, in um
    roiInfo = beaddata.roiInfo(iROI, 1:5);
    disp(['ROI# ' num2str(iROI) ' idx ' num2str(iTile) ' ij ' mat2str(ij)]);
    disp(mat2str(roiInfo));
    % = [centers radii ...
%                 (repmat(pxScale,nROIs,1) .* centers) + ...
%                 repmat(xyzThis(1:2)+imgOrig,nROIs,1)  ];
    dfi =  beadImages_now{ij(1),ij(2)}; %622 img = beadImages_now{ij(1),ij(2)}.img;
    axes(handles.axes2); cla; imshow(dfi(chan(1)).CDataS(greyLUT(1),greyLUT(2))); %622 imshow(uint8(255*mat2gray(img)));
    set(handles.axes2,'YDir',bsGUI.YDir2); % for Thorimage
    handles.axes2.YTick=[];
    if handles.cbShowID.Value
        viscircles(roiInfo(1:2),roiInfo(3),'LineWidth',1,'EdgeColor','m');
    end
    xyCenter = roiInfo(4:5);
    zOffset = (xyCenter-beadsurv.TL(1:2)) * beadsurv.dZdXY(:);
    xyzTarget = [xyCenter zOffset+beadsurv.TL(3)];
    disp(['Target is ' mat2str(xyzTarget)]);
    launchMove(xyzTarget);
    moveDone(xyzTarget);
    
end

function LUTArray = getLUTValues(handles)
LUTArray = [ str2double(handles.eTauLo.String) str2double(handles.eTauHi.String);
             0 str2double(handles.eLUTHi.String) ];

function [tZero,tZeroSaved] = processTZero(dfh,dfi,chans,handles)
uit = handles.uitTZero;

% do the fits
for k=1:numel(chans)
    ch = chans(k);
    [y,t] = dfh(ch).rangeData;
    fitter = DecayFitter_exp2_gauss;
    [vals,fail,tau,redchisq] = fitter.fitProgressive(t,y,0.07);
    p = fitter.params';
%     msg{1} = ['Fit with mean tau = ' num2str(tau,3) '  redChiSq = ' num2str(redchisq,3)];
%     msg{2} = mat2str(p,3);
    if isempty(redchisq)
        rChi(ch)  = fail;
    else
        rChi(ch)  = redchisq;
    end
    tZero(ch) = p(5);
    tZeroSaved(ch) = dfi(ch).tPeakIRF + dfi(k).tAdjust;
    uit.Data{k,1} = num2str(ch);
    uit.Data{k,2} = num2str(tZero(ch),3);
    if uit.Data{k,4}
        try 
            tZeroToFix = str2double(uit.Data{k,3});
            dfi(ch).setTPeak(tZeroToFix); 
        catch
            uit.Data{k,4} = false;
        end
    else
        uit.Data{k,3} = num2str(tZeroSaved(ch),3);
    end
end
% disp(' rChi tZero tZeroSaved');
% disp(mat2str([rChi(:) tZero(:) tZeroSaved(:)],3));

function uitTZero_ListenerCallback(src, eventdata, hObject, handles)
% hObject    handle to uitTZero (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
ij = eventdata.Indices;
if ij(2)==4 && hObject.Data{ij{1},4}
   % just set to true
   % respond to logical setting of fixed
   disp(hObject.Data{ij(1),ij(2)});
   hObject.Data{ij(1),ij(2)} = false;
end

% --- Executes when entered data in editable cell(s) in uitTZero.
function uitTZero_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitTZero (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
ij = eventdata.Indices;
if ij(2)==4
   if hObject.Data{ij(1),4}
       % just set to true
       % respond to logical setting of fixed
       hObject.Data{ij(1),3} = hObject.Data{ij(1),2};
   else
       hObject.Data{ij(1),3} = '';
   end
end 

% PRO-FORMA JUNK
function cbShowID_Callback(hObject, eventdata, handles)
function eRadiusShrink_Callback(hObject, eventdata, handles)
function eRadiusShrink_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eTauLo_Callback(hObject, eventdata, handles)
function eTauHi_Callback(hObject, eventdata, handles)
function cbLT_Callback(hObject, eventdata, handles)
function eLUTHi_Callback(hObject, eventdata, handles)
function eSmooth_Callback(hObject, eventdata, handles)
function eSmooth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function cbSmooth_Callback(hObject, eventdata, handles)
function cbRadiusShrink_Callback(hObject, eventdata, handles)
function eTestCol_Callback(hObject, eventdata, handles)
function eTestCol_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eTestRow_Callback(hObject, eventdata, handles)
function eTestRow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eTauLo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eTauHi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eLUTHi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function cbAutosaveBeaddata_Callback(hObject, eventdata, handles)
function eXYZ_BL_Callback(hObject, eventdata, handles)
function eXYZ_BR_Callback(hObject, eventdata, handles)
function eXYZ_TL_Callback(hObject, eventdata, handles)
function eXYZ_TR_Callback(hObject, eventdata, handles)
function eFileSpec_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eXYZ_BL_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eXYZ_BR_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eXYZ_TL_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eXYZ_TR_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eFullField_Callback(hObject, eventdata, handles)
% edit CreateFcn's
function eFFCenters_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eNumSubfields_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eSubfieldOffset_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eThreshold1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eBLCorner_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eSourceCh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eRadiusRange_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eSens_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eThreshold2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eSetNumber_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eFOV_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eFullField_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eSetNumber_Callback(hObject, eventdata, handles)
function eBLCorner_Callback(hObject, eventdata, handles)
function pbRedoFullField_Callback(hObject, eventdata, handles)
function eSourceCh_Callback(hObject, eventdata, handles)
function eRadiusRange_Callback(hObject, eventdata, handles)
function eThreshold1_Callback(hObject, eventdata, handles)
function eSens_Callback(hObject, eventdata, handles)
function eThreshold2_Callback(hObject, eventdata, handles)

