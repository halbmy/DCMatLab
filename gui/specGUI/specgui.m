function varargout = specgui(varargin)
% SPECGUI M-file for specgui.fig
%      SPECGUI, by itself, creates a new SPECGUI or raises the existing
%      singleton*.
%
%      H = SPECGUI returns the handle to a new SPECGUI or the handle to
%      the existing singleton*.
%
%      SPECGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPECGUI.M with the given input arguments.
%
%      SPECGUI('Property','Value',...) creates a new SPECGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before specgui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to specgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help specgui

% Last Modified by GUIDE v2.5 29-May-2008 14:01:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @specgui_OpeningFcn, ...
                   'gui_OutputFcn',  @specgui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before specgui is made visible.
function specgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to specgui (see VARARGIN)

% Choose default command line output for specgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes specgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = specgui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function disableall_Callback(hObject, eventdata, handles)
set(handles.checkbox1,'Enable','Off','Value',0);
set(handles.min1,'Enable','Off');
set(handles.max1,'Enable','Off');
set(handles.unit1,'Enable','Off');
set(handles.checkbox2,'Enable','Off','Value',0);
set(handles.min2,'Enable','Off');
set(handles.max2,'Enable','Off');
set(handles.unit2,'Enable','Off');
set(handles.checkbox3,'Enable','Off','Value',0);
set(handles.min3,'Enable','Off');
set(handles.max3,'Enable','Off');
set(handles.unit3,'Enable','Off');
set(handles.checkbox4,'Enable','Off','Value',0);
set(handles.min4,'Enable','Off');
set(handles.max4,'Enable','Off');
set(handles.unit4,'Enable','Off');
set(handles.checkbox5,'Enable','Off','Value',0);
set(handles.min5,'Enable','Off');
set(handles.max5,'Enable','Off');
set(handles.unit5,'Enable','Off');
set(handles.checkbox6,'Enable','Off','Value',0);
set(handles.min6,'Enable','Off');
set(handles.max6,'Enable','Off');
set(handles.unit6,'Enable','Off');
set(handles.checkbox7,'Enable','Off','Value',0);
set(handles.min7,'Enable','Off');
set(handles.max7,'Enable','Off');
set(handles.unit7,'Enable','Off');


% --- Executes on button press in loadfile.
function loadfile_Callback(hObject, eventdata, handles)
global A Unit Col filename
[fname,pname]=uigetfile({'*.asc';'*.*'},'Load data file');
if ~isstr(fname), return; end
specgui('disableall_Callback',gcbo,[],guidata(gcbo))       
filename=fullfile(pname,fname);
fid=fopen(filename);
Unit={};Col={};
zeile=destrip(fgetl(fid));
while isempty(zeile), zeile=destrip(fgetl(fid)); end
i=1;[Col{i},zeile]=strtok(zeile);
while ~isempty(zeile),
    i=i+1;
    [col,zeile]=strtok(zeile);
    if ~isempty(col), Col{i}=col; end
end
units=fgetl(fid);
i=1;[Unit{i},units]=strtok(units);
while ~isempty(units),
    i=i+1;
    [unit,units]=strtok(units);
    if ~isempty(unit), Unit{i}=unit; end
end
ss='';for i=1:length(Col), ss=[ss '%f']; end
A=fscanf(fid,'%f',[length(Col) Inf])';
% A=mytextscan(fid,ss);
fclose(fid);
A(A<=-999)=NaN;
set(handles.status,'String',sprintf('Read %s: %d data columns',fname,size(A,1)));
set(handles.zmin,'String',num2str(min(A(:,1))));
set(handles.zmax,'String',num2str(max(A(:,1))));
if size(A,2)>1,
    set(handles.checkbox1,'Enable','On','String',Col{2});
    set(handles.min1,'Enable','On','String',min(A(:,2)));
    set(handles.max1,'Enable','On','String',max(A(:,2)));
    set(handles.unit1,'Enable','On','String',Unit{2});
    set(handles.checkbox1,'Value',1);
end
if size(A,2)>2,
    set(handles.checkbox2,'Enable','On','String',Col{3});
    set(handles.min2,'Enable','On','String',min(A(:,3)));
    set(handles.max2,'Enable','On','String',max(A(:,3)));
    set(handles.unit2,'Enable','On','String',Unit{3});
    set(handles.checkbox2,'Value',1);
end
if size(A,2)>3,
    set(handles.checkbox3,'Enable','On','String',Col{4});
    set(handles.min3,'Enable','On','String',min(A(:,4)));
    set(handles.max3,'Enable','On','String',max(A(:,4)));
    set(handles.unit3,'Enable','On','String',Unit{4});
    set(handles.checkbox3,'Value',1);
end
if size(A,2)>4,
    set(handles.checkbox4,'Enable','On','String',Col{5});
    set(handles.min4,'Enable','On','String',min(A(:,5)));
    set(handles.max4,'Enable','On','String',max(A(:,5)));
    set(handles.unit4,'Enable','On','String',Unit{5});
    set(handles.checkbox4,'Value',1);
end
if size(A,2)>5,
    set(handles.checkbox5,'Enable','On','String',Col{6});
    set(handles.min5,'Enable','On','String',min(A(:,6)));
    set(handles.max5,'Enable','On','String',max(A(:,6)));
    set(handles.unit5,'Enable','On','String',Unit{6});
    set(handles.checkbox5,'Value',1);
end
if size(A,2)>6,
    set(handles.checkbox6,'Enable','On','String',Col{7});
    set(handles.min6,'Enable','On','String',min(A(:,7)));
    set(handles.max6,'Enable','On','String',max(A(:,7)));
    set(handles.unit6,'Enable','On','String',Unit{7});
    set(handles.checkbox6,'Value',1);
end
if size(A,2)>7,
    set(handles.checkbox7,'Enable','On','String',Col{8});
    set(handles.min7,'Enable','On','String',min(A(:,8)));
    set(handles.max7,'Enable','On','String',max(A(:,8)));
    set(handles.unit7,'Enable','On','String',Unit{8});
    set(handles.checkbox7,'Value',1);
end
set(handles.plotlogs,'Enable','On');
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       

% --- Executes during object creation, after setting all properties.
function min1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function min1_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes during object creation, after setting all properties.
function max1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function max1_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes during object creation, after setting all properties.
function min2_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function min2_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes during object creation, after setting all properties.
function max2_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function max2_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes during object creation, after setting all properties.
function min3_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function min3_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function edit6_Callback(hObject, eventdata, handles)


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes during object creation, after setting all properties.
function min4_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function min4_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       

% --- Executes during object creation, after setting all properties.
function max3_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function max3_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes during object creation, after setting all properties.
function max4_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function max4_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit8_Callback(hObject, eventdata, handles)


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes during object creation, after setting all properties.
function min5_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function min5_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes during object creation, after setting all properties.
function max5_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function max5_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       

% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes during object creation, after setting all properties.
function min6_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function min6_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       

% --- Executes during object creation, after setting all properties.
function max6_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function max6_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes during object creation, after setting all properties.
function min7_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function min7_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes during object creation, after setting all properties.
function max7_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function max7_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes during object creation, after setting all properties.
function zmin_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function zmin_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes during object creation, after setting all properties.
function zmax_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function zmax_Callback(hObject, eventdata, handles)
specgui('plotlogs_Callback',gcbo,[],guidata(gcbo))       


% --- Executes on button press in plotlogs.
function plotlogs_Callback(hObject, eventdata, handles)
global A B Unit Col cols
mindep=str2num(get(handles.zmin,'String'));
maxdep=str2num(get(handles.zmax,'String'));
cols=[];mincol=[];maxcol=[];
if get(handles.checkbox1,'Value'), 
    cols=[cols 1];
    mincol=[mincol str2num(get(handles.min1,'String'))];
    maxcol=[maxcol str2num(get(handles.max1,'String'))];
end
if get(handles.checkbox2,'Value'), 
    cols=[cols 2]; 
    mincol=[mincol str2num(get(handles.min2,'String'))];
    maxcol=[maxcol str2num(get(handles.max2,'String'))];
end
if get(handles.checkbox3,'Value'), 
    cols=[cols 3]; 
    mincol=[mincol str2num(get(handles.min3,'String'))];
    maxcol=[maxcol str2num(get(handles.max3,'String'))];
end
if get(handles.checkbox4,'Value'), 
    cols=[cols 4]; 
    mincol=[mincol str2num(get(handles.min4,'String'))];
    maxcol=[maxcol str2num(get(handles.max4,'String'))];
end
if get(handles.checkbox5,'Value'), 
    cols=[cols 5]; 
    mincol=[mincol str2num(get(handles.min5,'String'))];
    maxcol=[maxcol str2num(get(handles.max5,'String'))];
end
if get(handles.checkbox6,'Value'), 
    cols=[cols 6]; 
    mincol=[mincol str2num(get(handles.min6,'String'))];
    maxcol=[maxcol str2num(get(handles.max6,'String'))];
end
if get(handles.checkbox7,'Value'), 
    cols=[cols 7]; 
    mincol=[mincol str2num(get(handles.min7,'String'))];
    maxcol=[maxcol str2num(get(handles.max7,'String'))];
end
depth=A(:,1);
B=A((depth>=mindep)&(depth<maxdep),[1 cols+1]);
set(figure(1),'MenuBar','none','NumberTitle','off','Name','Log plot');clf;
ncols=length(cols);
colors={'b','r','g','c','m','y','k'};
for i=1:ncols,
    ncol=cols(i);
    logi=B(:,i+1);
    logi((logi<mincol(i))|(logi>maxcol(i)))=NaN;
    B(:,i+1)=logi;
    if ncols>1, subplot(1,ncols,i); end
    plot(logi,B(:,1),[colors{i} '-']);
    axis ij tight;grid on;
    xlabel([strrep(Col{ncol+1},'_',' ') ' in ' Unit{ncol+1}]);
    if i==1, ylabel('depth in m'); end
end
set(handles.makefft,'Enable','On');

% --- Executes on button press in makefft.
function makefft_Callback(hObject, eventdata, handles)
global B F
dx=abs(median(diff(B(:,1))));
F=0;%zeros(fix(size(B,1)/2),size(B,2));
for ncol=1:size(B,2)-1
    logi=B(:,ncol+1);
    fi1=find(isfinite(logi));    
    fi2=find(isnan(logi));
    logi(fi2)=interp1(B(fi1,1),logi(fi1),logi(fi2),'nearest');
    [rc,ic,ft]=fouriertransform(logi,dx);
    if ncol==1, F=ft(:); end    
    F(:,ncol+1)=sqrt(rc.^2+ic.^2);
end    
set(handles.fftaxis,'Enable','On');
set(handles.spectrum,'Enable','On');
set(handles.maxaxis,'Enable','On');
set(handles.ymin,'Enable','On');
set(handles.ymax,'Enable','On');

specgui('fftaxis_Callback',gcbo,[],guidata(gcbo))       


% --- Executes during object creation, after setting all properties.
function fftaxis_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in fftaxis.
function fftaxis_Callback(hObject, eventdata, handles)
global F
ty=get(handles.fftaxis,'Value');
if ty>1,
    set(handles.maxaxis,'String',num2str(1/F(2,1)));
else
    set(handles.maxaxis,'String',num2str(F(end,1)));    
end
specgui('spectrum_Callback',gcbo,[],guidata(gcbo))       


% --- Executes on button press in spectrum.
function spectrum_Callback(hObject, eventdata, handles)
global F cols Col
ncols=size(F,2)-1;
maxf=str2num(get(handles.maxaxis,'String'));
xla={'wavenumber in 1/m','wavelength in m'};
ty=get(handles.fftaxis,'Value');
figure(2);clf;
set(2,'NumberTitle','off','Name','Spectra');%,'MenuBar','none');
colors={'b','r','g','c','m','y','k'};
if ty>1, xma=min(find(F(:,1)>=1/maxf));
else xma=max(find(F(:,1)<=maxf)); end
for i=1:ncols,
    if ncols>1, subplot(ncols,1,i); end
    if ty>1,
        plot(1./(F(2:end,1)+eps),F(2:end,i+1),[colors{i} '.-']);
    else
        plot(F(2:end,1),F(2:end,i+1),[colors{i} '.-']);
    end
    set(gca,'xlim',[0 maxf]);
%     set(gca,'ylim',[0 max(F(2:xma,i+1))]);
%     set(gca,'ylim',[0 F(7,i+1)]);
    grid on;
    if i==ncols, xlabel(xla{ty}); end
    ylabel('amplitude');
    legend(strrep(Col{cols(i)+1},'_',' '));
%     legend(strrep(Col{ncol+1},'_',' '));
end
yl=get(gca,'YLim');
set(handles.ymin,'String',yl(1));
set(handles.ymax,'String',yl(2));



% --- Executes during object creation, after setting all properties.
function maxaxis_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function maxaxis_Callback(hObject, eventdata, handles)
specgui('spectrum_Callback',gcbo,[],guidata(gcbo))       


% --- Executes on button press in movingwindow.
function movingwindow_Callback(hObject, eventdata, handles)


% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
if strcmp(questdlg('Exit program?','Exit','Yes','No','Yes'),'Yes'),
    for i=20:-1:1, if ishandle(i), delete(i); end; end
    delete(gcbf);
end


% --- Executes on button press in pick.
function pick_Callback(hObject, eventdata, handles)
[x,y]=ginput(1)
hold on;plot(x,y,'ro');hold off;
uiwait(msgbox(sprintf('x=%g y=%g',x,y),'Position'));


% --- Executes on button press in zoom.
function zoom_Callback(hObject, eventdata, handles)
k=waitforbuttonpress;
cp=get(gca,'CurrentPoint'); 
finalRect=rbbox;
cp1=get(gca,'CurrentPoint'); 
set(gca,'Xlim',[max(min(cp(1,1),cp1(1,1)),0) max(cp(1,1),cp1(1,1))],...
    'Ylim',[min(cp(1,2),cp1(1,2)) max(cp(1,2),cp1(1,2))]);
% specgui('ymax_Callback',gcbo,[],guidata(gcbo))       


% --- Executes during object creation, after setting all properties.
function ymin_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ymin_Callback(hObject, eventdata, handles)
specgui('ymax_Callback',gcbo,[],guidata(gcbo))       


% --- Executes during object creation, after setting all properties.
function ymax_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ymax_Callback(hObject, eventdata, handles)
ymin=str2num(get(handles.ymin,'String'));
ymax=str2num(get(handles.ymax,'String'));
set(gca,'YLim',[ymin ymax]);


% --- Executes on button press in export.
function export_Callback(hObject, eventdata, handles)
global filename
[outfile,iseps,ispdf]=getgrfile(filename);
if iseps, % Testversion %n. Zeile einklammern
    epsprint(2,outfile,ispdf);
else 
    exportpng(2,outfile); 
end

