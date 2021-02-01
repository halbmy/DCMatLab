function varargout = wallbert(varargin)
% WALLBERT M-file for wallbert.fig
%      WALLBERT, by itself, creates a new WALLBERT or raises the existing
%      singleton*.
%
%      H = WALLBERT returns the handle to a new WALLBERT or the handle to
%      the existing singleton*.
%
%      WALLBERT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WALLBERT.M with the given input arguments.
%
%      WALLBERT('Property','Value',...) creates a new WALLBERT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before wallbert_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to wallbert_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help wallbert

% Last Modified by GUIDE v2.5 13-Mar-2008 11:25:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wallbert_OpeningFcn, ...
                   'gui_OutputFcn',  @wallbert_OutputFcn, ...
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

% --- Executes just before wallbert is made visible.
function wallbert_OpeningFcn(hObject, eventdata, handles, varargin)
global Meshopts Invopts
% Choose default command line output for wallbert
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using wallbert.
if exist('meshoptions.cfg','file'), 
    Meshopts=readconfigfile('meshoptions.cfg'); 
    if isfield(Meshopts,'meshtype')&&(Meshopts.meshtype>0), set(handles.meshtype,'Value',Meshopts.meshtype); end
    if isfield(Meshopts,'pararefine')&&(Meshopts.pararefine>0), set(handles.meshrefine,'Value',Meshopts.pararefine); end
end
if exist('invoptions.cfg','file'), Invopts=readconfigfile('invoptions.cfg'); end
if exist('wallbert.wbg','file'),
    AA=readconfigfile('wallbert.wbg');
    if isfield(AA,'geo'), set(handles.geometry,'Value',AA.geo); end
    numsides=[2 1 2 3 4];
    o2='Off';o3=o2;o4=o2;
    if numsides(AA.geo)>=2, o2='On'; end
    if numsides(AA.geo)>=3, o3='On'; end
    if numsides(AA.geo)>=4, o4='On'; end
    set(handles.n2edit,'Enable',o2);
    set(handles.d2edit,'Enable',o2);
    set(handles.n2text,'Enable',o2);
    set(handles.d2text,'Enable',o2);
    set(handles.n3edit,'Enable',o3);
    set(handles.d3edit,'Enable',o3);
    set(handles.n3text,'Enable',o3);
    set(handles.d3text,'Enable',o3);
    set(handles.n4edit,'Enable',o4);
    set(handles.d4edit,'Enable',o4);
    set(handles.n4text,'Enable',o4);
    set(handles.d4text,'Enable',o4);
    if isfield(AA,'l'), set(handles.ledit,'String',num2str(AA.l)); end
    if isfield(AA,'t'), set(handles.tedit,'String',num2str(AA.t)); end
    if isfield(AA,'a'), set(handles.aedit,'String',num2str(AA.a)); end
    if isfield(AA,'n1'), set(handles.n1edit,'String',num2str(AA.n1)); end
    if isfield(AA,'n2'), set(handles.n2edit,'String',num2str(AA.n2)); end
    if isfield(AA,'n3'), set(handles.n3edit,'String',num2str(AA.n3)); end
    if isfield(AA,'n4'), set(handles.n4edit,'String',num2str(AA.n4)); end
    if isfield(AA,'d1'), set(handles.d1edit,'String',num2str(AA.d1)); end
    if isfield(AA,'d2'), set(handles.d2edit,'String',num2str(AA.d2)); end
    if isfield(AA,'d3'), set(handles.d3edit,'String',num2str(AA.d3)); end
    if isfield(AA,'d4'), set(handles.d4edit,'String',num2str(AA.d4)); end
end
wallbert('draw_Callback',hObject,[],guidata(hObject));
% set(handles.loadgeometry,'Enable','on');
% set(handles.savegeometry,'Enable','on');

% UIWAIT makes wallbert wait for user response (see UIRESUME)
iconify(handles.figure1,'wallbert.ico');
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = wallbert_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function geometry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to geometry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in geometry.
function geometry_Callback(hObject, eventdata, handles)
numsides=[2 1 2 3 4];
geo=get(handles.geometry,'Value');
o2='Off';o3=o2;o4=o2;
if numsides(geo)>=2, o2='On'; end
if numsides(geo)>=3, o3='On'; end
if numsides(geo)>=4, o4='On'; end
set(handles.n2edit,'Enable',o2);
set(handles.d2edit,'Enable',o2);
set(handles.n2text,'Enable',o2);
set(handles.d2text,'Enable',o2);
set(handles.n3edit,'Enable',o3);
set(handles.d3edit,'Enable',o3);
set(handles.n3text,'Enable',o3);
set(handles.d3text,'Enable',o3);
set(handles.n4edit,'Enable',o4);
set(handles.d4edit,'Enable',o4);
set(handles.n4text,'Enable',o4);
set(handles.d4text,'Enable',o4);
wallbert('draw_Callback',gcbo,[],guidata(gcbo));


% --- Executes on button press in savegeometry.
function savegeometry_Callback(hObject, eventdata, handles)
global filename
AA.geo=get(handles.geometry,'Value');
AA.l=str2num(get(handles.ledit,'String'));
AA.t=str2num(get(handles.tedit,'String'));
AA.a=str2num(get(handles.aedit,'String'));
AA.n1=str2num(get(handles.n1edit,'String'));
AA.d1=str2num(get(handles.d1edit,'String'));
if AA.geo~=2, 
    AA.n2=str2num(get(handles.n2edit,'String'));
    AA.d2=str2num(get(handles.d2edit,'String'));
end
if AA.geo>3,
    AA.n3=str2num(get(handles.n3edit,'String'));
    AA.d3=str2num(get(handles.d3edit,'String'));
end
if AA.geo>4,
    AA.n4=str2num(get(handles.n4edit,'String'));
    AA.d4=str2num(get(handles.d4edit,'String'));
end
[pp,nn]=fileparts(filename);
infile=fullfile(pp,'*.wbg');
[fname,pname]=uiputfile('*.wbg','Save geometry file',infile);
[pp,nn,ee]=fileparts(fname);
if ~strcmp(lower(ee),'.wbg'), fname=[fname '.wgb']; end
if isstr(fname),
   writeconfigfile(fullfile(pname,fname),AA,'WallBERT Config File');
end

% --- Executes on button press in loadgeometry.
function loadgeometry_Callback(hObject, eventdata, handles)
global filename
[pp,nn]=fileparts(filename);
infile=fullfile(pp,'*.wbg');
[fname,pname]=uigetfile('*.wbg','Load geometry file',infile);
if ~isstr(fname), return; end
AA=readconfigfile(fullfile(pname,fname));
if isfield(AA,'geo'), set(handles.geometry,'Value',AA.geo); end
wallbert('geometry_Callback',gcbo,[],guidata(gcbo));
if isfield(AA,'l'), set(handles.ledit,'String',num2str(AA.l)); end
if isfield(AA,'t'), set(handles.tedit,'String',num2str(AA.t)); end
if isfield(AA,'a'), set(handles.aedit,'String',num2str(AA.a)); end
if isfield(AA,'n1'), set(handles.n1edit,'String',num2str(AA.n1)); end
if isfield(AA,'n2'), set(handles.n2edit,'String',num2str(AA.n2)); end
if isfield(AA,'n3'), set(handles.n3edit,'String',num2str(AA.n3)); end
if isfield(AA,'n4'), set(handles.n4edit,'String',num2str(AA.n4)); end
if isfield(AA,'d1'), set(handles.d1edit,'String',num2str(AA.d1)); end
if isfield(AA,'d2'), set(handles.d2edit,'String',num2str(AA.d2)); end
if isfield(AA,'d3'), set(handles.d3edit,'String',num2str(AA.d3)); end
if isfield(AA,'d4'), set(handles.d4edit,'String',num2str(AA.d4)); end
wallbert('draw_Callback',gcbo,[],guidata(gcbo));

% --- Executes during object creation, after setting all properties.
function ledit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function ledit_Callback(hObject, eventdata, handles)
wallbert('draw_Callback',gcbo,[],guidata(gcbo));


% --- Executes during object creation, after setting all properties.
function tedit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function tedit_Callback(hObject, eventdata, handles)
wallbert('draw_Callback',gcbo,[],guidata(gcbo));


% --- Executes during object creation, after setting all properties.
function aedit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function aedit_Callback(hObject, eventdata, handles)
wallbert('draw_Callback',gcbo,[],guidata(gcbo));


% --- Executes on button press in draw.
function draw_Callback(hObject, eventdata, handles)
geo=get(handles.geometry,'Value');
l=str2num(get(handles.ledit,'String'));
t=str2num(get(handles.tedit,'String'));
a=str2num(get(handles.aedit,'String'));
n1=str2num(get(handles.n1edit,'String'));
n2=str2num(get(handles.n2edit,'String'));
n3=str2num(get(handles.n3edit,'String'));
n4=str2num(get(handles.n4edit,'String'));
d1=str2num(get(handles.d1edit,'String'));
d2=str2num(get(handles.d2edit,'String'));
d3=str2num(get(handles.d3edit,'String'));
d4=str2num(get(handles.d4edit,'String'));
% ch=get(handles.axes1,'Children');    
% for i=1:length(ch), delete(ch(i)); end
axes(handles.axes1);cla;
line([0 1 1 0 0]*l,[0 0 1 1 0]*t,'Color','black');grid on;
% rectangle('Position',[0 0 l t]);
hold on;
d11=l-d1-(n1-1)*a;
green=[0 0.5 0];
%Box
text(0,t,'t','HorizontalAlignment','right','VerticalAlignment','middle','Color',green);
text(l,0,'l','VerticalAlignment','top','HorizontalAlignment','center','Color',green);
% 1st side always
x1=d11+((1:n1)-1)*a;
valid=(min(x1)>=0)&(max(x1)<=l);
plot(x1,0,'bx');
% for i=1:n1, plot(x1,0,'bx'); end
text(d11+0.5*a,0,'a','VerticalAlignment','bottom','HorizontalAlignment','center','Color',green);
text(l-d1/2,0,'d1','VerticalAlignment','bottom','HorizontalAlignment','center','Color','red');
text(d11,0,'1','VerticalAlignment','bottom','HorizontalAlignment','center','Color','blue');
if n1>1, text(d11+a,0,'2','VerticalAlignment','bottom','HorizontalAlignment','center','Color','blue'); end
text(l-d1,0,num2str(n1),'VerticalAlignment','bottom','HorizontalAlignment','center','Color','blue');
nsum=n1;
if geo==1, % opp. sides
    x2=l-d2-((1:n2)-1)*a;
    valid=valid&(min(x2)>=0)&(max(x2)<=l);
    plot(x2,t,'bx');
    text(l-d2/2,t,'d2','VerticalAlignment','top','HorizontalAlignment','center','Color','red');
    text(l-d2,t,num2str(n1+1),'VerticalAlignment','top','HorizontalAlignment','center','Color','blue');
    text(l-d2-(n2-1)*a,t,num2str(n1+n2),'VerticalAlignment','top','HorizontalAlignment','center','Color','blue');
    nsum=nsum+n2;
end
if geo>2, % 2,3 or 4 sides: side 2
    y2=d2+((1:n2)-1)*a;
    valid=valid&(min(y2)>=0)&(max(y2)<=t);
    plot(l,y2,'bx');
    text(l,d2/2,'d2','VerticalAlignment','middle','HorizontalAlignment','left','Color','red');
    text(l,d2,num2str(n1+1),'VerticalAlignment','middle','HorizontalAlignment','right','Color','blue');
    text(l,d2+(n2-1)*a,num2str(n1+n2),'VerticalAlignment','middle','HorizontalAlignment','right','Color','blue');
    nsum=nsum+n2;
end
if geo>3, % 3 or 4 sides: side 3
    x3=l-d3-((1:n3)-1)*a;
    valid=valid&(min(x3)>=0)&(max(x3)<=l);
    plot(x3,t,'bx');
    text(l-d3/2,t,'d3','VerticalAlignment','bottom','HorizontalAlignment','center','Color','red');
    text(l-d2,t,num2str(n1+n2+1),'VerticalAlignment','top','HorizontalAlignment','center','Color','blue');
    text(l-d2-(n2-1)*a,t,num2str(n1+n2+n3),'VerticalAlignment','top','HorizontalAlignment','center','Color','blue');
    nsum=nsum+n3;
end
if geo>4, % 4 sides: side 4
    y4=t-d4-((1:n4)-1)*a;
    valid=valid&(min(y4)>=0)&(max(y4)<=t);
    plot(0,y4,'bx');
    text(0,t-d4/2,0,'d4','VerticalAlignment','middle','HorizontalAlignment','left','Color','red');
    text(0,t-d4,num2str(n1+n2+n3+1),'VerticalAlignment','middle','HorizontalAlignment','left','Color','blue');
    text(0,t-d4-(n4-1)*a,num2str(n1+n2+n3+n4),'VerticalAlignment','middle','HorizontalAlignment','left','Color','blue');
    nsum=nsum+n4;
end
% set(gca,'XLim',[0 l],'Ylim',[0 t]);
hold off;
axis equal tight xy
sonoff='On';
if valid, 
    set(handles.status,'ForegroundColor','black','String',...
        sprintf('Set up %gx%gm geometry using %d electrodes',l,t,nsum));
else    
    sonoff='Off';
    set(handles.status,'ForegroundColor','red','String','Geometry invalid (electrodes outside)');
end
set(handles.loaddata,'Enable',sonoff);

    
% --- Executes during object creation, after setting all properties.
function n1edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function n1edit_Callback(hObject, eventdata, handles)
wallbert('draw_Callback',gcbo,[],guidata(gcbo));


% --- Executes during object creation, after setting all properties.
function n2edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function n2edit_Callback(hObject, eventdata, handles)
wallbert('draw_Callback',gcbo,[],guidata(gcbo));


% --- Executes during object creation, after setting all properties.
function n3edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function n3edit_Callback(hObject, eventdata, handles)
wallbert('draw_Callback',gcbo,[],guidata(gcbo));


% --- Executes during object creation, after setting all properties.
function n4edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function n4edit_Callback(hObject, eventdata, handles)
wallbert('draw_Callback',gcbo,[],guidata(gcbo));


% --- Executes during object creation, after setting all properties.
function d1edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function d1edit_Callback(hObject, eventdata, handles)
wallbert('draw_Callback',gcbo,[],guidata(gcbo));


% --- Executes during object creation, after setting all properties.
function d2edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function d2edit_Callback(hObject, eventdata, handles)
wallbert('draw_Callback',gcbo,[],guidata(gcbo));


% --- Executes during object creation, after setting all properties.
function d3edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function d3edit_Callback(hObject, eventdata, handles)
wallbert('draw_Callback',gcbo,[],guidata(gcbo));


% --- Executes during object creation, after setting all properties.
function d4edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function d4edit_Callback(hObject, eventdata, handles)
wallbert('draw_Callback',gcbo,[],guidata(gcbo));


%%%%%%%%%%%%%%%%%%%%DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in loaddata.
function loaddata_Callback(hObject, eventdata, handles)
global Data filename
alltypes='*.tx0;*.dat;*.TX0;*.DAT;*.ohm';
if isempty(filename), infile=alltypes; else
    [pname,fname,ext]=fileparts(filename);
    infile=strrep(filename,[fname ext],alltypes);
end 
[fname,pname]=uigetfile(alltypes,'Read data file',infile);
if ~isstr(fname), return; end
filename=fullfile(pname,fname);
if strcmp(filename(end-2:end),'tx0'), % Lippmann device file
    Data=readtx0file(filename);
    if isfield(Data,'sip'), 
        for i=1:length(Data.sip),
            sip=Data.sip{i};
            if isfield(sip,'i'), sip=delmeasurement(sip,sip.i==0); end
            if isfield(sip,'k'), sip=delmeasurement(sip,sip.k==0); end
            Data.sip{i}=sip;
        end
    end
    fi=find((Data.i==0)|(Data.u==0));
%     fi=find((Data.err>0.2)|(Data.i==0)|(abs(Data.u)<1e-3));
    Data=delmeasurement(Data,fi);
%     Data.allrho(fi,:)=[];
%     Data.allip(fi,:)=[];
    set(handles.status,'ForegroundColor','black','String',sprintf('Read TX0 file. %d data using %d electrodes',length(Data.a),size(Data.elec,1)));
else % other
    ch=check2dfile(filename);
    if ch==1, % res2dinvfile with wrong electrode positions (0,1,2...)
        Data=readres2dinvfile(filename);
        Data.rho=Data.r./Data.k;
        set(handles.status,'ForegroundColor','black','String',...
            sprintf('Read RES2DINV file. %d data using %d electrodes',length(Data.a),size(Data.elec,1)));
    elseif ch==2, % unified data format, see resistivity.net
        Data=readunifile(filename);
        if isfield(Data,'f'),
            Data.ff=unique(Data.f);
            if length(Data.ff)>1, % multiple frequencies
                for i=1:length(Data.ff),
                    fi=find(Data.f==Data.ff(i));
                    Data.sip{i}=rmfield(getmeasurement(Data,fi),'f');
                    Data.sip{i}.nr=fi;
                end
                fn=fieldnames(Data.sip{1});
                for i=1:length(fn), 
                    Data=setfield(Data,fn{i},getfield(Data.sip{1},fn{i})); 
                end
            end
        end
        if ~isfield(Data,'rho'),
            if isfield(Data,'i')&&isfield(Data,'u'), Data.rho=Data.u./Data.i;
            elseif isfield(Data,'r'),%
                kk=getkonf2d(Data);
                Data.rho=Data.r./kk;
            else
                errordlg('Could not get resistances from file',filename);return;
            end
        end
        set(handles.status,'ForegroundColor','black','String',...
            sprintf('Read unified format file. %d data using %d electrodes',length(Data.a),size(Data.elec,1)));        
    end
end
geo=get(handles.geometry,'Value');
n1=str2num(get(handles.n1edit,'String'));
n2=str2num(get(handles.n2edit,'String'));
n3=str2num(get(handles.n3edit,'String'));
n4=str2num(get(handles.n4edit,'String'));
nsum=n1+n2*(geo~=2)+n3*(geo>3)+n4*(geo>4);
if nsum~=size(Data.elec,1),
    set(handles.status,'ForegroundColor','red','String',...
        sprintf('Electrode size mismatch: %d (geometry) vs. %d (data)',nsum,size(Data.elec,1)));
    set(handles.makemesh,'Enable','Off');
    set(handles.calcgeo,'Enable','Off');
    set(handles.showdata,'Enable','Off');
    set(handles.inversion,'Enable','Off');
    set(handles.showmodel,'Enable','Off');
    set(handles.showphases,'Enable','Off');    
    return;
end
set(handles.figure1,'Name',['WALLBERT - ' strrep(filename,[pwd filesep],'')]);
l=str2num(get(handles.ledit,'String'));
t=str2num(get(handles.tedit,'String'));
a=str2num(get(handles.aedit,'String'));
d1=str2num(get(handles.d1edit,'String'));
d2=str2num(get(handles.d2edit,'String'));
d3=str2num(get(handles.d3edit,'String'));
d4=str2num(get(handles.d4edit,'String'));
Data.elec(1:n1,1)=l-d1-(n1-1:-1:0)*a;
Data.elec(1:n1,2)=0;
nsum=n1;
if geo==1,
    Data.elec(n1+1:n1+n2,1)=l-d2-(0:n2-1)*a;
    Data.elec(n1+1:n1+n2,2)=t;
    nsum=nsum+n2;
end
if geo>=3,
    Data.elec(n1+1:n1+n2,1)=l;
    Data.elec(n1+1:n1+n2,2)=d2+(0:n2-1)*a;
    nsum=nsum+n2;
end
if geo>=4,
    Data.elec(n1+n2+1:n1+n2+n3,1)=l-d3-(0:n3-1)*a;
    Data.elec(n1+n2+1:n1+n2+n3,2)=t;
    nsum=nsum+n3;
end
if geo>=5,
    Data.elec(n1+n2+n3+1:n1+n2+n3+n4,1)=0;
    Data.elec(n1+n2+n3+1:n1+n2+n3+n4,2)=t-d2-(0:n4-1)*a;
    nsum=nsum+n4;
end
if size(Data.elec,2)>2, Data.elec(:,3:end)=[]; end
% DD=Data;
% if isfield(DD,'r'), DD=rmfield(DD,'r'); end
% if isfield(DD,'err'), DD=rmfield(DD,'err'); end
% saveinv2dfile('tomo.ohm',DD,1,'#x y');    
% % if isfield(Data,'u'), saveinv2dfile('tomo.ohm',rmfield(Data,{'r','rho','err'}),1,'#x y');
% % else saveinv2dfile('tomo.ohm',rmfield(Data,{'r'}),1,'#x y'); end
Data.x=Data.elec(:,1);Data.y=Data.elec(:,2);
Data.elec(1:nsum,1)=1:nsum;Data.elec(:,2)=0;
set(handles.makemesh,'Enable','On');
set(handles.calcgeo,'Enable','On');
set(handles.showdata,'Enable','Off');
set(handles.inversion,'Enable','Off');
set(handles.showmodel,'Enable','Off');
set(handles.showphases,'Enable','Off');
wallbert('initfields_Callback',gcbo,[],guidata(gcbo))

% --- Executes on button press in calcgeo.
function calcgeo_Callback(hObject, eventdata, handles)
global Mesh Data Meshopts
set(handles.status,'ForegroundColor','black','String','Making mesh for primary potential...');
set(handles.figure1,'Pointer','Watch');pause(0.1);
% Da=readunifile('tomo.ohm');
% eigenes Netz
l=str2num(get(handles.ledit,'String'));
t=str2num(get(handles.tedit,'String'));
a=str2num(get(handles.aedit,'String'));
Poly=[];
Poly.edge=[1 2;2 3;3 4;4 1];Poly.edge(:,3)=-1;
Poly.node=[0 0;l 0;l t;0 t];
primrefine=0.1;primquality=33.3;primmaxcellsize=0;
if isfield(Meshopts,'primlocalrefine')&&~isempty(Meshopts.primlocalrefine)&&(Meshopts.primlocalrefine>0), primrefine=Meshopts.primlocalrefine; end
if isfield(Meshopts,'primquality')&&~isempty(Meshopts.primquality)&&(Meshopts.primquality>=25)&&(Meshopts.primquality<35), primquality=Meshopts.primquality; end
if isfield(Meshopts,'primmaxcellsize')&&~isempty(Meshopts.primmaxcellsize)&&(Meshopts.primmaxcellsize>0), primmaxcellsize=Meshopts.primmaxcellsize; end
dr=a*primrefine;
el1=[Data.x Data.y];
mid=mean(el1);    
for i=1:size(el1,1), 
    del=mid-el1(i,:);
    del=del/norm(del);
    el1(i,:)=el1(i,:)+del*dr;
end
Poly.node=[Poly.node;Data.x Data.y];
Poly.node(:,3)=-99;Poly.node(1:4,3)=0;
el1(:,3)=0;
Poly.node=[Poly.node;el1];
if 1,
    Poly.node(end+1,:)=[l*0.25 t/2 -999];
    Poly.node(end+1,:)=[l*0.75 t/2 -1000];
else
    Poly.node(1,3)=-999;
    Poly.node(1,4)=-1000;
end
Poly.region=[mid 2 primmaxcellsize];
writepoly2d('tomoPrim.poly',Poly);
dos(['dctriangle -q' num2str(primquality) ' tomoPrim.poly']);
MeshPrim=loadmesh('tomoPrim.bms');
if ~isfield(MeshPrim,'nnodes')||(MeshPrim.nnodes<100),
    set(handles.status,'ForegroundColor','red','String','Meshing cancelled');
    errordlg('Something went wrong creating mesh','Meshing error');
    set(handles.figure1,'Pointer','Arrow');
    return;
end
set(handles.status,'ForegroundColor','black','String',...
    sprintf('Read primary Mesh (%d nodes). Calculating...',MeshPrim.nnodes));
pause(0.2);
tripatchmod(MeshPrim);pause(0.5);
fid=fopen('rho1.map','w');
fprintf(fid,'0 1\n1 1\n2 1\n');
fclose(fid);
%old: H3 mesh
% fid=fopen('rho.map','w');
% for i=2:Mesh.ncells+1, fprintf(fid,'%d 1\n',i); end
% % fprintf(fid,'1 1\n2 1\n');
% fclose(fid);
% dos('meshconvert -vd2 -r2 -o tomoPrim -YB tomo.bms');

if 1, %dcfemlib2, nur collect-file
%     dos('dcmod -v -o num -a rho.map -s tomo.ohm tomoPrim.bms');
%     Num=readunifile('num.ohm');
%     dos('dcmod -v -o num -a rho1.map -p pot\prim\pot tomoPrim.bms');
    if 1,%dcmod
        dos('dcmod -v -o num -a rho1.map tomoPrim.bms');
        MEA=readcollect('num.collect',1);
    else
%         dos('dcfem -v -d2 -s2 -D -C -o num -a rho1.map tomoPrim.bms');
        dos('meshconvert -vd2 -D -p -o tomoPrim2.bms tomoPrim.bms');
        dos('dcfem -v -d2 -s2 -D -C -o num -a rho1.map tomoPrim2.bms');
        MEA=readcollect('num.collect');
    end
    Data.k=ones(size(Data.a));
    [u1,u2]=collectrhoa(Data,MEA);
    Data.k=1./u1;
    Data.r=abs(Data.rho.*Data.k);
    if isfield(Data,'sip'),
       for i=1:length(Data.sip),
           if isfield(Data.sip{i},'u')&isfield(Data.sip{i},'i'), 
               tmpdata=Data.sip{i};
               tmpdata.k=ones(size(tmpdata.a));
               [u1,u2]=collectrhoa(tmpdata,MEA);
               Data.sip{i}.k=1./u1;
               Data.sip{i}.r=tmpdata.u./tmpdata.i./u1;
           end
       end
    end
else, %dcfemlib1, provides pot/prim/pot.D.*->(interpolate)->pot/prim/
    dos('dcfem -v -C -D -d2 -s2 -o pot/prim/pot tomoPrim.bms');
    MEA=readcollect('pot.collect');
    Num=readunifile('tomo.ohm');
    Num.k=ones(size(Num.a));
    [u1,u2]=collectrhoa(Num,MEA);
    Num.k=1./u1;
    Num.r=Num.rho.*Num.k; % ??? useless....
end
% Data.k=1./Num.rho;
% Data.r=abs(Data.rho.*Data.k);
set(handles.figure1,'Pointer','Arrow');
if isfield(Data,'ip'),
    set(handles.status,'ForegroundColor','black','String',...
        sprintf('min/max rhoa = %g/%g Ohmm, min/max phase = %g/%g mrad',...
        rndig(min(Data.r)),rndig(max(Data.r)),rndig(min(Data.ip)),rndig(max(Data.ip))));
else
    set(handles.status,'ForegroundColor','black','String',...
        sprintf('min/max rhoa = %g/%g Ohmm',rndig(min(Data.r)),rndig(max(Data.r))));    
end
set(handles.showdata,'Enable','On');
wallbert('initfields_Callback',gcbo,[],guidata(gcbo))
wallbert('showdata_Callback',gcbo,[],guidata(gcbo))

% --- Executes on selection change in datafield.
function initfields_Callback(hObject, eventdata, handles)
global Data
if isfield(Data,'sip')&&isfield(Data,'ff'),
    set(handles.choosef,'Enable','On','String',num2strcell(Data.ff,'%gHz'));
end
ss={};
if isfield(Data,'r'), ss{end+1}='rhoa'; end
if isfield(Data,'ip'), ss{end+1}='IP'; end
if isfield(Data,'k'), ss{end+1}='k'; end
if isfield(Data,'i'), ss{end+1}='I'; end
if isfield(Data,'u'), ss{end+1}='U'; 
elseif isfield(Data,'rho'), ss{end+1}='R'; end
if isfield(Data,'err'), ss{end+1}='err'; end
set(handles.datafield,'String',ss,'Value',1);
set(handles.editmin,'Enable','on');
set(handles.mintext,'Enable','on');
set(handles.editmax,'Enable','on');
set(handles.maxtext,'Enable','on');
set(handles.trim,'Enable','on');

% --- Executes on selection change in choosef.
function choosef_Callback(hObject, eventdata, handles)
global Data
nfr=get(handles.choosef,'Value');
if isfield(Data,'sip')&&(length(Data.sip)>=nfr);
    ff=fieldnames(Data.sip{nfr});
    for i=1:length(ff),
        Data=setfield(Data,ff{i},getfield(Data.sip{nfr},ff{i}));
    end
end
wallbert('initfields_Callback',gcbo,[],guidata(gcbo))
wallbert('showdata_Callback',gcbo,[],guidata(gcbo))

% --- Executes during object creation, after setting all properties.
function datafield_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in datafield.
function datafield_Callback(hObject, eventdata, handles)
wallbert('showdata_Callback',gcbo,[],guidata(gcbo))

% --- Executes on button press in showdata.
function showdata_Callback(hObject, eventdata, handles)
global Data
axes(handles.axes1);cla;
nf=get(handles.datafield,'Value');
ss=get(handles.datafield,'String');
sakt=ss{nf};
% if strcmp(lower(get(handles.choosef,'Enable')),'on')&&isfield(Data,'sip'),%spectral data
%     nfr=get(handles.choosef,'Value');
%     if length(Data.sip)>=nfr, tmpdata=Data.sip{nfr}; end
% end

field=[];canot='';
if isfield(Data,'r'), field=Data.r;canot='Ohmm'; end
if isequal(lower(sakt),'k'), field=Data.k;canot='m'; end
if isequal(lower(sakt),'u'), field=abs(Data.u)*1000;canot='mV'; end
if isequal(lower(sakt),'r'), field=Data.rho;canot='Ohmm'; end
if isequal(lower(sakt),'i'), field=abs(Data.i)*1000;canot='mA'; end
if isequal(lower(sakt),'ip'), field=Data.ip;canot='mrad'; end
if isequal(lower(sakt),'err'), field=Data.err*100;canot='%'; end

if get(handles.showhist,'Value'),
    if min(field)>0, loghist(field); else hist(field,30); end
    xlabel(canot);ylabel('histogram');
else
    opt=struct('canot',canot,'cauto',1,'clog',min(field)>0);
    [Data.mids,Data.seps,Data.ii,Data.kk]=showdata2d(Data,field,opt);
    im=get(handles.axes1,'Children');
    set(im,'ButtonDownFcn','wallbert(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))');
end
set(handles.editmin,'String',num2str(rndig(min(field)*0.99)));
set(handles.editmax,'String',num2str(rndig(max(field)*1.01)));


% --- Executes on button press in showhist.
function showhist_Callback(hObject, eventdata, handles)
wallbert('showdata_Callback',gcbo,[],guidata(gcbo))


% --- Executes during object creation, after setting all properties.
function editmin_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editmin_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function editmax_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editmax_Callback(hObject, eventdata, handles)


% --- Executes on button press in trim.
function trim_Callback(hObject, eventdata, handles)
global Data
nf=get(handles.datafield,'Value');
ss=get(handles.datafield,'String');
sakt=ss{nf};
field=Data.r;
if isequal(lower(sakt),'k'), field=Data.k; end
if isequal(lower(sakt),'u'), field=abs(Data.u)*1000; end
if isequal(lower(sakt),'i'), field=abs(Data.i)*1000; end
if isequal(lower(sakt),'ip'), field=Data.ip; end

mi=str2num(get(handles.editmin,'String'));
ma=str2num(get(handles.editmax,'String'));
fi=find((field<mi)|(field>ma));
ndata=length(Data.r);
if isempty(fi),
    set(handles.status,'ForegroundColor','red','String','Nothing to delete!');
elseif length(fi)==ndata,
    set(handles.status,'ForegroundColor','red','String','Trying to delete all data! Cancelled...');
else
    Data=delmeasurement(Data,fi);
    set(handles.status,'ForegroundColor','black','String',...
        sprintf('Deleted %d values ...',length(fi)));
    if isfield(Data,'sip'),
        for i=1:length(Data.sip),
            sipdata=Data.sip{i};
            field=sipdata.r;
            if isfield(sipdata,'k')&&isequal(lower(sakt),'k'), field=sipdata.k; end
            if isfield(sipdata,'u')&&isequal(lower(sakt),'u'), field=abs(sipdata.u)*1000; end
            if isfield(sipdata,'u')&&isequal(lower(sakt),'i'), field=abs(sipdata.i)*1000; end
            if isfield(sipdata,'u')&&isequal(lower(sakt),'ip'), field=sipdata.ip; end
            fi=find((field<mi)|(field>ma));
            if ~isempty(fi),
                sipdata=delmeasurement(sipdata,fi);
                Data.sip{i}=sipdata;
            end
        end
    end
end
wallbert('showdata_Callback',gcbo,[],guidata(gcbo))


%%%%%%%%%%%%%%%%%%%%%%%%%MESH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function meshtype_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in meshtype.
function meshtype_Callback(hObject, eventdata, handles)
global Meshopts
Meshopts.meshtype=get(handles.meshtype,'Value');

% --- Executes on button press in meshrefine.
function meshrefine_Callback(hObject, eventdata, handles)
global Meshopts
Meshopts.pararefine=get(handles.meshrefine,'Value');

% --- Executes on button press in makemesh.
function makemesh_Callback(hObject, eventdata, handles)
global Data Mesh Sortdata Meshopts
% geo=get(handles.geometry,'Value');
l=str2num(get(handles.ledit,'String'));
t=str2num(get(handles.tedit,'String'));
a=str2num(get(handles.aedit,'String'));
xi=unique(Data.x);
xi(xi==0)=[];
xi(xi==l)=[];
x=[fillgap(xi(1),a) xi(:)' l-fliplr(fillgap(l-xi(end),a))];
yi=unique(Data.y);
if length(yi)>2,
    yi(yi==0)=[];
    yi(yi==t)=[];
    y=[fillgap(yi(1),a) yi(:)' t-fliplr(fillgap(l-yi(end),a))];
else
    y=linspace(0,t,round(t/a)+1);
end
% save xy x y
type=get(handles.meshtype,'Value')-1;
refine=get(handles.meshrefine,'Value');
switch type,
    case 1, % irregularly spaced rectangular mesh        
        if refine,
            xx=[x;x+[diff(x)/2 0]];
            x=xx(1:end-1);
            yy=[y;y+[diff(y)/2 0]];
            y=yy(1:end-1);
        end
        Mesh=grid4mesh(x,y);
    case 2,% regularly spaced rectangular mesh
        x=linspace(0,l,round(l/a*(1+refine))+1);
        y=linspace(0,t,round(t/a*(1+refine))+1);
        Mesh=grid4mesh(x,y);
    case 3, % unstructured triangular mesh
        Poly=[];Poly.nnodes=length(Data.x)+6;
        Poly.node=[0 0;l 0;l t;0 t;Data.x Data.y;l*0.25 t/2;l*0.75 t/2];
        Poly.nodemarker=ones(Poly.nnodes,1)*(-99);
        Poly.nodemarker(1:4)=0;Poly.nodemarker(end-1:end)=-1000:-999;
        Poly.node(:,3)=Poly.nodemarker;
        Poly.edge=[1 2;2 3;3 4;4 1];Poly.edge(:,3)=-1;
        paramaxcellsize=0;
        if isfield(Meshopts,'paramaxcellsize')&&(Meshopts.paramaxcellsize>=0),
            paramaxcellsize=Mesh.paramaxcellsize;
        elseif refine, paramaxcellsize=a*a/2.5; end
        Poly.region=[l/2 t/2 2 paramaxcellsize];
        writepoly2d('tomo.poly',Poly);
        paraquality=33.3;
        if isfield(Meshopts,'paraquality')&&(Meshopts.paraquality>=25), paraquality=Meshopts.paraquality; end
        dos(sprintf('dctriangle -v -q%g tomo.poly',paraquality));
        Mesh=loadmesh('tomo.bms');
    otherwise, % rectangular mesh refined to triangles
        if refine,
            xx=[x;x+[diff(x)/2 0]];
            x=xx(1:end-1);
            yy=[y;y+[diff(y)/2 0]];
            y=yy(1:end-1);
        end
        Mesh=grid2mesh(ones(length(x)-1,length(y)-1),x,y);
end
Mesh.node=abs(Mesh.node);
axes(handles.axes1);cla;
tripatchmod(Mesh);
hold on;plot(Data.x,Data.y,'ro');hold off;% axis xy;

% fi=find(Mesh.nodemarker==-99);
% hold on;plot(Mesh.node(fi,1),Mesh.node(fi,2),'ro');hold off
if type>0,
    fi999=find((Mesh.node(:,1)==0)&(Mesh.node(:,2)==0));
    fi1000=find((Mesh.node(:,1)==max(Mesh.node(:,1)))&(Mesh.node(:,2)==0));
else % secmesh refine overwrites 999/1000 nodemarker when on boundary
    fi999=111;while Mesh.nodemarker(fi999)<0, fi999=fi999+1; end
    fi1000=222;while Mesh.nodemarker(fi1000)<0, fi1000=fi1000+1; end
end
Mesh.nodemarker(fi999)=-999;
Mesh.nodemarker(fi1000)=-1000;
Mesh.cellattr(:)=(2:Mesh.ncells+1);
savemesh(Mesh,'tomo.bms');
if 0, %dcfemlib2->interpolate
    dos('interpolate -v -d2 -C -D -i pot/prim/pot -q pot/sec/pot tomoPrim.bms meshSec.bms');
end
set(handles.inversion,'Enable','On');
set(handles.status,'String',sprintf('Created mesh with %d cells and %d nodes',Mesh.ncells,Mesh.nnodes));

% --- Executes on button press in inversion.
function inversion_Callback(hObject, eventdata, handles)
global Data Mesh Model Invopts
type=get(handles.meshtype,'Value')-1;
if type<2,
    Sortdata=rmfield(Data,{'r','k'});
    for i=1:length(Data.x),
        fff=find((Data.x(i)==Mesh.node(:,1))&(Data.y(i)==Mesh.node(:,2)));
        if ~isempty(fff), 
            fin(i)=fff;
            Mesh.nodemarker(fff)=-99; 
            dosave=1;
        end
    end
    savemesh(Mesh,'tomo.bms');
    [so,ind]=sort(fin);
    Sortdata.elec=[Data.x(ind) Data.y(ind)];
    aind=1:length(ind);aind(ind)=1:length(ind);
    Sortdata.a=aind(Data.a);Sortdata.b=aind(Data.b);
    Sortdata.m=aind(Data.m);Sortdata.n=aind(Data.n);
    saveinv2dfile('tomo.ui',Sortdata,1,'#x y');
else
    saveinv2dfile('tomo.ui',Data,1,'#x y');
end
set(handles.status,'ForegroundColor','black','String','Running inversion...');
set(handles.figure1,'Pointer','Watch');pause(0.1);
cmdline='dcinv -v -p tomo.bms';
if isfield(Invopts,'withprimpot')&&(Invopts.withprimpot>0), cmdline=[cmdline ' -x pot/sec/pot']; end
if isfield(Invopts,'robustdata')&&(Invopts.robustdata>0), cmdline=[cmdline ' -R']; end
if isfield(Invopts,'blockymodel')&&(Invopts.blockymodel>0), cmdline=[cmdline ' -B']; end
if isfield(Invopts,'optlambda')&&(Invopts.optlambda>0), cmdline=[cmdline ' -O -l 1000']; 
elseif isfield(Invopts,'lambda'), cmdline=[cmdline ' -l ' num2str(Invopts.lambda)]; end
if isfield(Invopts,'lowerbound')&&(Invopts.lowerbound>0), cmdline=[cmdline ' -b ' num2str(Invopts.lowerbound)]; end
if isfield(Invopts,'upperbound')&&(Invopts.upperbound>0), cmdline=[cmdline ' -u ' num2str(Invopts.upperbound)]; end
if isfield(Invopts,'maxiter')&&(Invopts.maxiter>0), cmdline=[cmdline ' -i ' num2str(Invopts.maxiter)]; end
cmdline=[cmdline ' tomo.ui'];
dos(cmdline);
if exist('resistivity.vec','file'),
    fid=fopen('resistivity.vec');res=fscanf(fid,'%f');fclose(fid);
    if length(res)==Mesh.ncells,
        Model.res=res;
        set(handles.showmodel,'Enable','On');
        wallbert('showmodel_Callback',gcbo,[],guidata(gcbo))       
    else
        errordlg('Size mismatches');
    end
else
    errordlg('File resistivity.vec not found','File not found');
end
if exist('phase.vec','file'),
    fid=fopen('phase.vec');phase=fscanf(fid,'%f');fclose(fid);
    if length(phase)==Mesh.ncells, 
        set(handles.showphases,'Enable','On'); 
        Model.phase=phase;
    end
end
set(handles.figure1,'Pointer','Arrow');pause(0.1);
set(handles.status,'ForegroundColor','black','String',...
    sprintf('min/max rho = %g/%g Ohmm, min/max phase = %g/%g mrad',...
    rndig(min(res)),rndig(max(res)),rndig(min(phase)),rndig(max(phase))));
if isfield(Data,'sip'),
    wallbert('sipinversion_Callback',hObject,[],guidata(hObject));
end

% --- Executes on button press in showmodel.
function showmodel_Callback(hObject, eventdata, handles)
global Data Mesh Model
if isfield(Model,'res')&isfield(Mesh,'ncells'),
    if length(Model.res)==Mesh.ncells,
        axes(handles.axes1);cla;
        dx=max(Mesh.node(:,1))-min(Mesh.node(:,1));
        dy=max(Mesh.node(:,2))-min(Mesh.node(:,2));
        [cmin,cmax]=tripatchmod(Mesh,Model.res,struct('cauto',1,'canot','Ohmm','cbar',1+(dx<2*dy)));
    else
        errordlg('Size mismatches!');
    end
end
% set(handles.phasnr,'Enable','On','String',num2strcell(Data.ff,'%gHz'));
im=get(handles.axes1,'Children');
set(im,'ButtonDownFcn','wallbert(''axes1_ButtonDownModel'',gcbo,[],guidata(gcbo))');


% --- Executes on button press in showphases.
function showphases_Callback(hObject, eventdata, handles)
global Data Mesh Model
if isfield(Model,'phase'),
    axes(handles.axes1);cla;
    dx=max(Mesh.node(:,1))-min(Mesh.node(:,1));
    dy=max(Mesh.node(:,2))-min(Mesh.node(:,2));
    [cmin,cmax]=tripatchmod(Mesh,Model.phase,struct('cauto',1,'canot','mrad','cbar',1+(dx<2*dy)));
end
im=get(handles.axes1,'Children');
set(im,'ButtonDownFcn','wallbert(''axes1_ButtonDownModel'',gcbo,[],guidata(gcbo))');


% ---------------------------MOUSE-BUTTON-----------------------------------
function axes1_ButtonDownFcn(hObject, eventdata, handles) 
global Data
cp=get(handles.axes1,'CurrentPoint'); 
nk=round(cp(1,2));
[aa,nx]=min(abs(Data.mids-cp(1,1)));
nr=find((Data.kk==nk)&(Data.ii==nx));
if ~isempty(nr),
    if length(nr)>1, nr=nr(end); end
    mess=sprintf('Nr %d: rhoa=%.1f',nr,Data.r(nr));
    if isfield(Data,'ip')&&(~isempty(Data.ip))&&(length(Data.ip)>=nr), 
        mess=sprintf('%s IP=%.2fmrad',mess,Data.ip(nr)); end  
    if isfield(Data,'k')&&(length(Data.k)>=nr), 
        mess=sprintf('%s k=%.2gm',mess,Data.k(nr)); end  
    if isfield(Data,'i')&&(length(Data.i)>=nr), 
        mess=sprintf('%s I=%.1fmA',mess,Data.i(nr)*1000); end  
    if isfield(Data,'u')&&(length(Data.u)>=nr), 
        mess=sprintf('%s U=%.1fmV',mess,Data.u(nr)*1000); end  
    if isfield(Data,'err')&&(~isempty(Data.err))&&(length(Data.err)>=nr), 
        mess=sprintf('%s Error=%.1f%%',mess,Data.err(nr)*100); end  
    mess=sprintf('%s\na=%d (%.1fm)',mess,Data.a(nr),Data.elec(Data.a(nr),1));
    if Data.b(nr)>0, mess=sprintf('%s  b=%d (%.1fm)',mess,Data.b(nr),Data.elec(Data.b(nr),1)); end
    mess=sprintf('%s  m=%d (%.1fm)',mess,Data.m(nr),Data.elec(Data.m(nr),1));
    if Data.n(nr)>0, mess=sprintf('%s  n=%d (%.1fm)',mess,Data.n(nr),Data.elec(Data.n(nr),1)); end
    if strcmp(questdlg(mess,'Delete Datum?','Yes','No','No'),'Yes'),
        Data=delmeasurement(Data,nr);
        sss=sprintf('min/max rhoa = %g/%g Ohmm',rndig(min(Data.r)),rndig(max(Data.r)));
        if isfield(Data,'ip')&&~isempty(Data.ip),
            sss=sprintf('%s, min/max phase = %g/%g mrad',rndig(min(Data.ip)),rndig(max(Data.ip)));
        end
        set(handles.status,'ForegroundColor','black','String',sss);        
        wallbert('showdata_Callback',gcbo,[],guidata(gcbo))       
    end
end


% --- Executes during object creation, after setting all properties.
function phasnr_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in phasnr.
function phasnr_Callback(hObject, eventdata, handles)
global Model
nr=get(handles.phasnr,'Value');
if isfield(Model,'resistivities')&&iscell(Model.resistivities)&&(length(Model.resistivities)>=nr),
    Model.res=Model.resistivities{nr}; end
if isfield(Model,'phases')&&iscell(Model.phases)&&(length(Model.phases)>=nr),
    Model.phase=Model.phases{nr}; end
wallbert('showmodel_Callback',gcbo,[],guidata(gcbo))       
% wallbert('showphases_Callback',gcbo,[],guidata(gcbo))       


%%%%%%%%%%%%%%%%%%%MENU%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No'), return; end
delete(handles.figure1)

% --------------------------------------------------------------------
function HelpMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function AboutMenu_Callback(hObject, eventdata, handles)
hh=msgbox({'WallBERT';'Boundless Electrical Resistivity Tomography';'on Walls';...
        'version 0.9';'';'Authos: T. Günther (& C. Rücker)';'resistivity.net production'},'About');
iconify(hh,'wallbert.ico');
uiwait(hh);


% --------------------------------------------------------------------
function releaseMenu_Callback(hObject, eventdata, handles)
ibrowse(['file://' pwd filesep '%23release.txt']);

% --------------------------------------------------------------------
function LoadGeometryMenuItem_Callback(hObject, eventdata, handles)
wallbert('loadgeometry_Callback',gcbo,[],guidata(gcbo))       

% --------------------------------------------------------------------
function LoadDataMenuItem_Callback(hObject, eventdata, handles)
wallbert('loaddata_Callback',gcbo,[],guidata(gcbo))       

% --------------------------------------------------------------------
function SaveDataMenuItem_Callback(hObject, eventdata, handles)
global Data filename
[pp,nn,ee]=fileparts(filename);
infile=fullfile(pp,[nn '.dat']);
[fname,pname]=uiputfile('*.dat ','Save data file',infile);
if isstr(fname),
    Savedata=Data;
    if isfield(Data,'sip'), % put together all the frequencies
        Savedata=combinesipdata(Data);
    end
    [pp,nn,ee]=fileparts(fname);
    if isempty(ee), fname=[fname '.dat']; end
    filename=fullfile(pname,fname);
    saveunifile(filename,Savedata);%,1,'#x y');
end

% --------------------------------------------------------------------
function SaveGeometryMenuItem_Callback(hObject, eventdata, handles)
wallbert('savegeometry_Callback',gcbo,[],guidata(gcbo))       


% --------------------------------------------------------------------
function OptionsMenu_Callback(hObject, eventdata, handles)
set(handles.GraphicsOptionsMenuItem,'Enable','off');
% set(handles.InversionOptionsMenuItem,'Enable','off');
% set(handles.MeshOptionsMenuItem,'Enable','off');

% --------------------------------------------------------------------
function MeshOptionsMenuItem_Callback(hObject, eventdata, handles)
global Meshopts
uiwait(meshoptions);
if isfield(Meshopts,'meshtype'), set(handles.meshtype,'Value',Meshopts.meshtype); end
if isfield(Meshopts,'pararefine'), set(handles.meshrefine,'Value',Meshopts.pararefine); end

% --------------------------------------------------------------------
function InversionOptionsMenuItem_Callback(hObject, eventdata, handles)
uiwait(invoptions);

% --------------------------------------------------------------------
function DefaultGeometryMenuItem_Callback(hObject, eventdata, handles)
AA.geo=get(handles.geometry,'Value');
AA.l=str2num(get(handles.ledit,'String'));
AA.t=str2num(get(handles.tedit,'String'));
AA.a=str2num(get(handles.aedit,'String'));
AA.n1=str2num(get(handles.n1edit,'String'));
AA.d1=str2num(get(handles.d1edit,'String'));
if AA.geo~=2, 
    AA.n2=str2num(get(handles.n2edit,'String'));
    AA.d2=str2num(get(handles.d2edit,'String'));
end
if AA.geo>3,
    n3=str2num(get(handles.n3edit,'String'));
    d3=str2num(get(handles.d3edit,'String'));
end
if AA.geo>4,
    n4=str2num(get(handles.n4edit,'String'));
    d4=str2num(get(handles.d4edit,'String'));
end
writeconfigfile('wallbert.wbg',AA,'WallBERT Config File');


% --- Executes during object creation, after setting all properties.
function choosef_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --------------------------------------------------------------------
function exportmenu_Callback(hObject, eventdata, handles)
global Model Mesh
smod=sonoff(isfield(Model,'res')&&isfield(Mesh,'ncells')&&(length(Model.res)==Mesh.ncells));
set(handles.exportgraphicsMenu,'Enable',smod);
set(handles.exportmodelMenu,'Enable',smod);
set(handles.exportspectrumMenu,'Enable',sonoff(ishandle(10)));

% --------------------------------------------------------------------
function exportgraphicsMenu_Callback(hObject, eventdata, handles)
global filename Model Mesh
f=figure(9);set(9,'Menubar','none');
% set(f,'Color',[1 1 1]);clf;
% pcolor(rand(3));
% cb=colorbar('horiz');
% set(gca,'Units','normalized');
% set(cb,'Units','normalized');
% set(handles.axes1,'Units','pixel');
% poall=get(handles.figure1,'Position');
% copyobj(handles.axes1,9);
mal=struct('cauto',1,'canot','Ohmm');%,'cbar',1+(dx<2*dy));
tripatchmod(Mesh,Model.res,mal);
if isempty(filename), filename='new.dat'; end
[pp,nn,ee]=fileparts(filename);
[fname,pname]=uiputfile('*.eps','Export Figure',strrep(filename,ee,'.eps'));
if ischar(fname),
    epsprint(9,fullfile(pname,fname),1);
end
if ishandle(9), delete(9); end

% --------------------------------------------------------------------
function exportmodelMenu_Callback(hObject, eventdata, handles)
global Mesh Data Model filename
% ascii file with midpoint
A=zeros(Mesh.ncells,2);
for i=1:Mesh.ncells, A(i,:)=mean(Mesh.node(Mesh.cell(i,:),:)); end
firstline='#x\ty';formstr='%.3f\t%.3f';
if isfield(Model,'resistivities')&&(length(Model.resistivities)>1),
    nres=length(Model.resistivities);ff=1:nres;
    if isfield(Data,'ff')&&(length(Data.ff)==nres), ff=Data.ff; end
    for i=1:nres,
        if length(Model.resistivities{i})==Mesh.ncells,
            A(:,end+1)=Model.resistivities{i};
            firstline=sprintf('%s\tres(%gHz)',firstline,ff(i));
            formstr=[formstr '\t%.4g'];
        end
        if isfield(Model,'phases')&&(length(Model.phases)>=i)&&(length(Model.phases{i})==Mesh.ncells),
            A(:,end+1)=Model.phases{i};
            firstline=sprintf('%s\tphi(%gHz)',firstline,ff(i));
            formstr=[formstr '\t%.4g'];
        end
    end
else
    if isfield(Model,'res')&&(length(Model.res)==Mesh.ncells), 
        A(:,end+1)=Model.res(:);firstline=[firstline '\tres/Ohmm']; 
        formstr=[formstr '\t%.4g']; end
    if isfield(Model,'phase')&&(length(Model.phase)==Mesh.ncells), 
        A(:,end+1)=Model.phase(:);firstline=[firstline '\tphase/mrad']; 
        formstr=[formstr '\t%.4g']; end
end
firstline=[firstline '\r\n'];formstr=[formstr '\r\n'];
[pp,nn,ee]=fileparts(filename);
infile=fullfile(pp,[nn '.out']);
[fname,pname]=uiputfile({'*.out';'*.mod';'*.dat';'*.*'},'Save output file',infile);
if isstr(fname),
    outname=fullfile(pname,fname);
    fid=fopen(outname,'w');
    if fid<0, errordlg('Could not open file!');return; end
    fprintf(fid,firstline);
    fprintf(fid,formstr,A');
    fclose(fid);
    set(handles.status,'ForegroundColor','black','String',...
        sprintf('Exported %d values to %s',size(A,1),outname));
end

% --------------------------------------------------------------------
function pickcellMenu_Callback(hObject, eventdata, handles)
im=get(handles.axes1,'Children');
set(im,'ButtonDownFcn','wallbert(''axes1_ButtonDownModel'',gcbo,[],guidata(gcbo))');


% --------------------------------------------------------------------
function sipinversion_Callback(hObject, eventdata, handles)
global Mesh Model Data Invopts
S=loadprimpot('Sens.mat');
abmn=[Data.a(:) Data.b(:) Data.m(:) Data.n(:)];
umin=50e-6;lambda=30;errperc=3;
if isfield(Data,'u'), errperc=1; end
if isfield(Invopts,'lambda')&&(Invopts.lambda>0), lambda=Invopts.lambda; end
if isfield(Invopts,'errperc')&&(Invopts.errperc>0), errperc=Invopts.errperc; end
if isfield(Invopts,'umin')&&(Invopts.umin>0), umin=Invopts.umin/1000; end

if isfield(Data,'u'), err=errperc+umin./Data.u; else err=errperc; end
C=onesidesmoothness(Mesh);CTC=C'*C;
xzero=zeros(Mesh.ncells,1);
for i=1:length(Data.sip),
    data=Data.sip{i}; % do data selection
    abmni=[data.a(:) data.b(:) data.m(:) data.n(:)];
    if size(abmn,2)~=size(abmni,2),
        fprintf('column size mismatch: %d vs. $d\n',size(abmn,2),size(abmni,2));
    end
    [tf,loc]=ismember(abmni,abmn,'Rows');
%     isequal(abmn(loc(tf),:),abmni(tf,:))
    dR=log(abs(data.r(tf)))-log(abs(Data.r(loc(tf))));
    aa=1./log(err(loc(tf))+1);
    aa(data.r(tf)<0)=0;
    aa(Data.r(loc(tf))<0)=0;
    D=spdiags(aa(:),0,sum(tf),sum(tf));
    dm=cglscdp(S(loc(tf),:),dR,lambda,CTC,D,1,xzero,xzero,100,1);
%     DSi=D*S(loc(tf),:);
%     dm=(DSi'*DSi+lam*CTC)\((D*dR(:))'*DSi)';
%     Model.phases{i}=(DSi'*DSi+lam*CTC)\((D*ipdata(:))'*DSi)';
    Model.resistivities{i} = Model.res.*exp(dm);
    ipdata=data.ip(tf);
    Model.phases{i}=cglscdp(S(loc(tf),:),ipdata(:),lambda,CTC,D,1,xzero,xzero,100,1);
end
set(handles.phasnr,'Enable','On','String',num2strcell(Data.ff,'%gHz'));
wallbert('phasnr_Callback',hObject,[],guidata(hObject));


function axes1_ButtonDownModel(hObject, eventdata, handles) 
global Mesh Model Data
cp=get(handles.axes1,'CurrentPoint');
if ~isfield(Mesh,'cellmids')||(size(Mesh.cellmids,1)~=Mesh.ncells),
    Mesh.cellmids=zeros(Mesh.ncells,Mesh.dim);
    for i=1:Mesh.ncells, Mesh.cellmids(i,:)=mean(Mesh.node(Mesh.cell(i,:),:)); end
end
di=sqrt((Mesh.cellmids(:,1)-cp(1,1)).^2+(Mesh.cellmids(:,2)-cp(1,2)).^2);
[mindi,nc]=min(di);
outstr=sprintf('x=%.2fm,y=%.2fm',Mesh.cellmids(nc,:));
if isfield(Model,'res')&&(length(Model.res)>=nc),
    outstr=sprintf('%s,res=%.2f Ohmm',outstr,Model.res(nc)); end
if isfield(Model,'phase')&&(length(Model.phase)>=nc),
    outstr=sprintf('%s,phase=%.2fmrad',outstr,Model.phase(nc)); end
if isfield(Model,'resistivities')&&iscell(Model.resistivities),
    nres=length(Model.resistivities);fax=1:nres;
    if isfield(Data,'ff')&&(length(Data.ff)==nres), fax=Data.ff; end
    fname=sprintf('Cell %d: %s',nc,outstr);
    set(figure(10),'MenuBar','none','NumberTitle','off','Name',fname);clf;
    isphases=isfield(Model,'phases')&&(length(Model.phases)==nres);
    if isphases, subplot(2,1,1); end
    resis=zeros(nres,1);
    for i=1:nres, resis(i)=Model.resistivities{i}(nc); end
    semilogx(fax,resis,'x-');title('amplitudes');grid on;
    set(gca,'XTick',fax,'XTickLabel',num2strcell(fax));
    xlabel('f in Hz');ylabel('res. in Ohmm');
    if isphases,
        subplot(2,1,2);phasis=zeros(nres,1);
        for i=1:nres, phasis(i)=Model.phases{i}(nc); end
        semilogx(fax,phasis,'x-');title('phases');grid on;
        set(gca,'XTick',fax,'XTickLabel',num2strcell(fax));
        xlabel('f in Hz');ylabel('phase in mrad');
    end
else
    uiwait(msgbox(outstr,['Picked Cell: ' num2str(nc)])); 
end


% --------------------------------------------------------------------
function exportspectrumMenu_Callback(hObject, eventdata, handles)
global filename
if ishandle(10),
    if isempty(filename), filename='new.dat'; end
    [pp,nn,ee]=fileparts(filename);
    [fname,pname]=uiputfile('*.eps','Export Figure',strrep(filename,ee,'-spec.eps'));
    if ischar(fname),
        epsprint(handles.figure1,fullfile(pname,fname),1);
    end
end

function erg=sonoff(value)
erg='Off';
if (nargin>0)&&isnumeric(value)&&(value>0), erg='On'; end
if (nargin>0)&&islogical(value)&&(value), erg='On'; end 


% --------------------------------------------------------------------
function exportfigureMenu_Callback(hObject, eventdata, handles)
global filename
if isempty(filename), filename='new.dat'; end
[pp,nn,ee]=fileparts(filename);
[fname,pname]=uiputfile('*.eps','Export Figure',strrep(filename,ee,'-fig.eps'));
if ischar(fname),
    epsprint(handles.figure1,fullfile(pname,fname),1);
end


% --------------------------------------------------------------------
function SaveOptionsMenuItem_Callback(hObject, eventdata, handles)
global Meshopts Invopts
writeconfigfile('meshoptions.cfg',Meshopts,'Mesh options');
writeconfigfile('invoptions.cfg',Invopts,'Inversion options');

% --------------------------------------------------------------------
function GraphicsOptionsMenuItem_Callback(hObject, eventdata, handles)


