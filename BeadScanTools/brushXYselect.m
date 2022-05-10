function varargout = brushXYselect(varargin)
% BRUSHXYSELECT MATLAB code for brushXYselect.fig
%      BRUSHXYSELECT, by itself, creates a new BRUSHXYSELECT or raises the existing
%      singleton*.
%
%      H = BRUSHXYSELECT returns the handle to a new BRUSHXYSELECT or the handle to
%      the existing singleton*.
%
%      BRUSHXYSELECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BRUSHXYSELECT.M with the given input arguments.
%
%      BRUSHXYSELECT('Property','Value',...) creates a new BRUSHXYSELECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before brushXYselect_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to brushXYselect_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help brushXYselect

% Last Modified by GUIDE v2.5 23-Feb-2021 08:19:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @brushXYselect_OpeningFcn, ...
                   'gui_OutputFcn',  @brushXYselect_OutputFcn, ...
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


% --- Executes just before brushXYselect is made visible.
function brushXYselect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to brushXYselect (see VARARGIN)

% Choose default command line output for brushXYselect
handles.output = hObject;
global beaddata
scatter(handles.axes1, ...
    beaddata.roiInfo(:,4),beaddata.roiInfo(:,5),1);
handles.brush = brush(handles.figure1);
handles.brush.ActionPostCallback = @(f,a) brushingDone(f,a);
handles.brushedPositions = [];
handles.figure1.UserData = handles;
brush(handles.figure1,'on');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes brushXYselect wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = brushXYselect_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.brushedPositions;
close(handles.figure1);

function brushingDone(src,eventdata)
global beaddata
ax = eventdata.Axes;
plt = ax.Children(1);
bdata = plt.BrushData;
src.UserData.brushedPositions = find(bdata);
if isempty(src.UserData.brushedPositions)
    str = 'No points selected';
else
    str = [num2str(sum(bdata)) ' points selected'];
end
src.UserData.tMessage.String = str;

function pbSelect_Callback(hObject, eventdata, handles)
handles.brushedPositions = ...
    handles.figure1.UserData.brushedPositions;
guidata(hObject,handles);
uiresume(handles.figure1);

function pbCancel_Callback(hObject, eventdata, handles)
handles.brushedPositions = [];
guidata(hObject,handles);
uiresume(handles.figure1);

