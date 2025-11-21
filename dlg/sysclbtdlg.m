function varargout = sysclbtdlg(varargin)
% SYSCLBTDLG MATLAB code for sysclbtdlg.fig
%      SYSCLBTDLG, by itself, creates a new SYSCLBTDLG or raises the existing
%      singleton*.
%
%      H = SYSCLBTDLG returns the handle to a new SYSCLBTDLG or the handle to
%      the existing singleton*.
%
%      SYSCLBTDLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SYSCLBTDLG.M with the given input arguments.
%
%      SYSCLBTDLG('Property','Value',...) creates a new SYSCLBTDLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sysclbtdlg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sysclbtdlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sysclbtdlg

% Last Modified by GUIDE v2.5 03-Oct-2025 18:36:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sysclbtdlg_OpeningFcn, ...
                   'gui_OutputFcn',  @sysclbtdlg_OutputFcn, ...
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


% --- Executes just before sysclbtdlg is made visible.
function sysclbtdlg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sysclbtdlg (see VARARGIN)

% Choose default command line output for sysclbtdlg
handles.output = hObject;
movegui(gcf,'center');
% Update handles structure
guidata(hObject, handles);
init_gui(hObject, handles);

% UIWAIT makes sysclbtdlg wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function init_gui(hObject, handles)
glvs
if ~isfield(handles,'setting')
    % lon0 lat0 hgt0 g0 inerN
    handles.file = 'c:\ygm\imu.txt';
    handles.setting = [34; 106; 400; 9.78; 5];
    handles.clbtresstr = num2str([eye(3); zeros(1,3); eye(3); zeros(1,3)]);
end
setting = handles.setting;
set(handles.simufile, 'String', handles.file);
set(handles.lon0, 'String', setting(1));
set(handles.lat0, 'String', setting(2));
set(handles.hgt0, 'String', setting(3));
set(handles.g0, 'String', setting(4));
set(handles.iterN, 'String', setting(5));
set(handles.clbtres, 'String', handles.clbtresstr);
guidata(hObject, handles);

function setting_change(handles, hObject, value, idx)
    if isnan(value), uiwait(msgbox('请输入一数值！','modal')); set(hObject,'String',handles.setting(idx)); return; end
    if value<0, uiwait(msgbox('输入数值不能为负！','modal')); set(hObject,'String',handles.setting(idx)); return; end
    if value>1200 && idx==8, uiwait(msgbox('总飞行时间须小于1200s！','modal')); set(hObject,'String',handles.setting(idx)); return; end
    handles.setting(idx) = value;
    guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = sysclbtdlg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function simufile_Callback(hObject, eventdata, handles)
% hObject    handle to simufile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of simufile as text
%        str2double(get(hObject,'String')) returns contents of simufile as a double
    setting_change(handles, hObject, str2double(get(hObject,'String')), 1);


% --- Executes during object creation, after setting all properties.
function simufile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simufile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function lon0_Callback(hObject, eventdata, handles)
% hObject    handle to lon0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lon0 as text
%        str2double(get(hObject,'String')) returns contents of lon0 as a double
    setting_change(handles, hObject, str2double(get(hObject,'String')), 1);


% --- Executes during object creation, after setting all properties.
function lon0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lon0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function lat0_Callback(hObject, eventdata, handles)
% hObject    handle to lat0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
    setting_change(handles, hObject, str2double(get(hObject,'String')), 2);


% --- Executes during object creation, after setting all properties.
function lat0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lat0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function hgt0_Callback(hObject, eventdata, handles)
% hObject    handle to hgt0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hgt0 as text
%        str2double(get(hObject,'String')) returns contents of hgt0 as a double
    setting_change(handles, hObject, str2double(get(hObject,'String')), 3);


% --- Executes during object creation, after setting all properties.
function hgt0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hgt0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function g0_Callback(hObject, eventdata, handles)
% hObject    handle to g0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of g0 as text
%        str2double(get(hObject,'String')) returns contents of g0 as a double
    setting_change(handles, hObject, str2double(get(hObject,'String')), 4);


% --- Executes during object creation, after setting all properties.
function g0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to g0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function iterN_Callback(hObject, eventdata, handles)
% hObject    handle to iterN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iterN as text
%        str2double(get(hObject,'String')) returns contents of iterN as a double
    setting_change(handles, hObject, str2double(get(hObject,'String')), 5);


% --- Executes during object creation, after setting all properties.
function iterN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iterN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in selectfile.
function selectfile_Callback(hObject, eventdata, handles)
% hObject    handle to simulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( {'*.txt'}, '选择惯组文件');
if isequal(filename,0), return; end
cd(pathname);
handles.file = [pathname,filename];
init_gui(hObject, handles);
guidata(hObject, handles);

% --- Executes on button press in simulate.
function simulate_Callback(hObject, eventdata, handles)
% hObject    handle to simulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dph = pi/180/3600;  ug = 9.8/1e6;
imu = load(handles.file);
imuplot(imu);
pos0 = posset(handles.setting(1), handles.setting(2), handles.setting(3));
clbt = sysclbt(imu, pos0, handles.setting(4), eye(3), handles.setting(5));
handles.clbtresstr = num2str([clbt.Kg; clbt.eb'/dph; clbt.Ka; clbt.db'/ug]);
init_gui(hObject, handles);
guidata(hObject, handles);


function clbtres_Callback(hObject, eventdata, handles)
% hObject    handle to clbtres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of clbtres as text
%        str2double(get(hObject,'String')) returns contents of clbtres as a double


% --- Executes during object creation, after setting all properties.
function clbtres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clbtres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
