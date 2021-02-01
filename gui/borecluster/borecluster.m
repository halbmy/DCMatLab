function varargout = borecluster(varargin)
% BORECLUSTER M-file for borecluster.fig
%      BORECLUSTER, by itself, creates a new BORECLUSTER or raises the existing
%      singleton*.
%
%      H = BORECLUSTER returns the handle to a new BORECLUSTER or the handle to
%      the existing singleton*.
%
%      BORECLUSTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BORECLUSTER.M with the given input arguments.
%
%      BORECLUSTER('Property','Value',...) creates a new BORECLUSTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before borecluster_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to borecluster_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help borecluster

% Last Modified by GUIDE v2.5 11-Oct-2007 09:42:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @borecluster_OpeningFcn, ...
                   'gui_OutputFcn',  @borecluster_OutputFcn, ...
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


% --- Executes just before borecluster is made visible.
function borecluster_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to borecluster (see VARARGIN)

% Choose default command line output for borecluster
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes borecluster wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = borecluster_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadfile.
function loadfile_Callback(hObject, eventdata, handles)
global A dep filename col
if isempty(filename), filename='*.asc'; end
[fname,pname]=uigetfile('*.asc;*.ASC','Read ASCII Data File',filename);
if ~isstr(fname), return; end
filename=fullfile(pname,fname);
hl=0;
fid=fopen(filename);
while isempty(destrip(fgetl(fid))), hl=hl+1; end
fclose(fid);
A=textread(filename,'','headerlines',hl+2);
A(A<-990)=NaN;
dep=A(:,1);
dep(end+1)=2*dep(end)-dep(end-1);
%%
fid=fopen(filename,'r');
for i=1:hl+1, zeile=fgetl(fid); end
fclose(fid);
i=1;col={};
while ~isempty(zeile), 
    [coli,zeile]=strtok(zeile);
    if isempty(coli), break; end
    col{i}=coli;
    i=i+1; 
end
nc=length(col);
%%
figure(11);na=size(A,2)-1;
for i=1:na,
    subplot(1,na,i);
    fi=find(A(:,i+1)>-900);
    plot(A(fi,i+1),dep(fi));
    set(gca,'YDir','reverse');
    if i<length(col), title(col{i+1}); end
end
%%
for i=1:size(A,2),
    aa=A(:,i);ma=max(aa);mi=min(aa);
    fi0=find(isnan(aa));
    fi1=find(~isnan(aa));
    aa(fi0)=interp1(dep(fi1),aa(fi1),dep(fi0),'linear','extrap');
    if mi==ma, aa(:)=0; else aa=(aa-mi)/(ma-mi); end
    A(:,i)=aa;
end
if nc>=1, set(handles.checkbox1,'String',col{1},'Enable','On','Value',0); 
else set(handles.checkbox1,'Enable','Off'); end
if nc>=2, set(handles.checkbox2,'String',col{2},'Enable','On','Value',1); 
else set(handles.checkbox2,'Enable','Off'); end
if nc>=3, set(handles.checkbox3,'String',col{3},'Enable','On','Value',1); 
else set(handles.checkbox3,'Enable','Off'); end
if nc>=4, set(handles.checkbox4,'String',col{4},'Enable','On','Value',1); 
else set(handles.checkbox4,'Enable','Off'); end
if nc>=5, set(handles.checkbox5,'String',col{5},'Enable','On','Value',1); 
else set(handles.checkbox5,'Enable','Off'); end
if nc>=6, set(handles.checkbox6,'String',col{6},'Enable','On','Value',1); 
else set(handles.checkbox6,'Enable','Off'); end
if nc>=7, set(handles.checkbox7,'String',col{7},'Enable','On','Value',1); 
else set(handles.checkbox7,'Enable','Off'); end
if nc>=8, set(handles.checkbox8,'String',col{8},'Enable','On','Value',1); 
else set(handles.checkbox8,'Enable','Off'); end
if nc>=9, set(handles.checkbox9,'String',col{9},'Enable','On','Value',1); 
else set(handles.checkbox9,'Enable','Off'); end
if nc>=10, set(handles.checkbox10,'String',col{10},'Enable','On','Value',1); 
else set(handles.checkbox10,'Enable','Off'); end
if nc>=11, set(handles.checkbox11,'String',col{11},'Enable','On','Value',1); 
else set(handles.checkbox11,'Enable','Off'); end
if nc>=12, set(handles.checkbox12,'String',col{12},'Enable','On','Value',1); 
else set(handles.checkbox12,'Enable','Off'); end
if nc>=13, set(handles.checkbox13,'String',col{13},'Enable','On','Value',1); 
else set(handles.checkbox13,'Enable','Off'); end
if nc>=14, set(handles.checkbox14,'String',col{14},'Enable','On','Value',1); 
else set(handles.checkbox14,'Enable','Off'); end
if nc>=15, set(handles.checkbox15,'String',col{15},'Enable','On','Value',1); 
else set(handles.checkbox15,'Enable','Off'); end
if nc>=16, set(handles.checkbox16,'String',col{16},'Enable','On','Value',1); 
else set(handles.checkbox16,'Enable','Off'); end
if nc>=17, set(handles.checkbox17,'String',col{17},'Enable','On','Value',1); 
else set(handles.checkbox17,'Enable','Off'); end
if nc>=18, set(handles.checkbox18,'String',col{18},'Enable','On','Value',1); 
else set(handles.checkbox18,'Enable','Off'); end
if nc>=19, set(handles.checkbox19,'String',col{19},'Enable','On','Value',1); 
else set(handles.checkbox19,'Enable','Off'); end
if nc>=20, set(handles.checkbox20,'String',col{20},'Enable','On','Value',1); 
else set(handles.checkbox20,'Enable','Off'); end
if nc>=21, set(handles.checkbox21,'String',col{21},'Enable','On','Value',1); 
else set(handles.checkbox21,'Enable','Off'); end
if nc>=22, set(handles.checkbox22,'String',col{22},'Enable','On','Value',1); 
else set(handles.checkbox22,'Enable','Off'); end
if nc>=23, set(handles.checkbox23,'String',col{23},'Enable','On','Value',1); 
else set(handles.checkbox23,'Enable','Off'); end
if nc>=24, set(handles.checkbox24,'String',col{24},'Enable','On','Value',1); 
else set(handles.checkbox24,'Enable','Off'); end
if nc>=25, set(handles.checkbox25,'String',col{25},'Enable','On','Value',1); 
else set(handles.checkbox25,'Enable','Off'); end
if nc>=26, set(handles.checkbox26,'String',col{26},'Enable','On','Value',1); 
else set(handles.checkbox26,'Enable','Off'); end
if nc>=27, set(handles.checkbox27,'String',col{27},'Enable','On','Value',1); 
else set(handles.checkbox27,'Enable','Off'); end
if nc>=28, set(handles.checkbox28,'String',col{28},'Enable','On','Value',1); 
else set(handles.checkbox28,'Enable','Off'); end
if nc>=29, set(handles.checkbox29,'String',col{29},'Enable','On','Value',1); 
else set(handles.checkbox29,'Enable','Off'); end
if nc>=30, set(handles.checkbox30,'String',col{30},'Enable','On','Value',1); 
else set(handles.checkbox30,'Enable','Off'); end
set(handles.cluster,'Enable','off');
set(handles.export,'Enable','off');



% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11


% --- Executes on button press in checkbox12.
function checkbox12_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox12


% --- Executes on button press in checkbox13.
function checkbox13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox13


% --- Executes on button press in checkbox14.
function checkbox14_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox14


% --- Executes on button press in checkbox15.
function checkbox15_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox15


% --- Executes on button press in checkbox16.
function checkbox16_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox16


% --- Executes on button press in checkbox17.
function checkbox17_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox17


% --- Executes on button press in checkbox18.
function checkbox18_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox18


% --- Executes on button press in checkbox19.
function checkbox19_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox19


% --- Executes on button press in checkbox20.
function checkbox20_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox20


% --- Executes on button press in checkbox21.
function checkbox21_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox21


% --- Executes on button press in checkbox22.
function checkbox22_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox22


% --- Executes on button press in checkbox23.
function checkbox23_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox23


% --- Executes on button press in checkbox24.
function checkbox24_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox24


% --- Executes on button press in checkbox25.
function checkbox25_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox25




% --- Executes on button press in checkbox26.
function checkbox26_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox26


% --- Executes on button press in checkbox27.
function checkbox27_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox27


% --- Executes on button press in checkbox28.
function checkbox28_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox28


% --- Executes on button press in checkbox29.
function checkbox29_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox29


% --- Executes on button press in checkbox30.
function checkbox30_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox30


% --- Executes on button press in dist.
function dist_Callback(hObject, eventdata, handles)
global A LI
set(handles.figure1,'Pointer','watch'); 
spalt=[];i=0;
i=i+1;if get(handles.checkbox1,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox2,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox3,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox4,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox5,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox6,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox7,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox8,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox9,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox10,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox11,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox12,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox13,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox14,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox15,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox16,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox17,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox18,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox19,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox20,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox21,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox22,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox23,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox24,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox25,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox26,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox27,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox28,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox29,'Value'), spalt(end+1)=i; end
i=i+1;if get(handles.checkbox30,'Value'), spalt(end+1)=i; end
% for i=1:30,
%     ii=eval(['get(handles.checkbox' num2str(i) ',''Value'')']);
%     if ii>0, spalt=[spalt i]; end
% end
meth='seuclidean';
% ss=get(handles.distmethod,'String');
% nn=get(handles.distmethod,'Value');
% meth=ss{nn}
DI=pdist(A(:,spalt),meth);
meth='ward';
% ss=get(handles.linkmethod,'String');
% nn=get(handles.linkmethod,'Value');
% meth=ss{nn}
LI=linkage(DI,meth);
set(handles.cluster,'Enable','on');
set(handles.figure1,'Pointer','arrow'); 
figure(1);clf;set(1,'Color',[1 1 1]);
dendrogram(LI);

function nclust_Callback(hObject, eventdata, handles)
% hObject    handle to nclust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nclust as text
%        str2double(get(hObject,'String')) returns contents of nclust as a double


% --- Executes during object creation, after setting all properties.
function nclust_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nclust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in cluster.
function cluster_Callback(hObject, eventdata, handles)
global LI dep filename
cl=str2num(get(handles.nclust,'String'));
CL=cluster(LI,'maxclust',cl);
%%
zm = dep(1:end-1)+diff(dep)/2;
A=[zm CL];
[pp,nn,ext]=fileparts(filename);
fid=fopen(strrep(filename,ext,['-cl' num2str(cl)]));
fprintf(fid,'%f %d\n',A');
fclose(fid);
%xlswrite(strrep(filename,ext,['-cl' num2str(cl)]),A);
% xlswrite(strrep(filename,ext,['-cl' num2str(cl) '.xls']),A);
%%
%colors={'r','g','b','c','m','y','k'};
colors={'r','g','b','c','m','y',[.5 .5 .5],[1 .5 .5],[.5 1 .5],[.5 .5 1],[.5 1 1],[1 .5 1],[1 1 .5]};
figure(1);clf;set(1,'Color',[1 1 1]);
for i=1:length(CL),
    col=colors{mod(CL(i)-1,length(colors))+1};
%     if CL(i)==7, col=[1 1 1]/2; end
    patch([-1 1 1 -1]*5,[dep(i) dep(i) dep(i+1) dep(i+1)],col,'EdgeColor',col);
end
axis equal tight
box on
set(gca,'YDir','reverse','XTick',[]);
yl=get(gca,'Ylim');dy=50;
for i=round(yl(1)/dy)*dy:10:round(yl(2)/dy)*dy, line([-5 -3],[i i],'Color','black'); end
set(handles.export,'Enable','on');


% --- Executes on button press in export.
function export_Callback(hObject, eventdata, handles)
global filename
[pp,nn,ee]=fileparts(filename);
if ishandle(1),
    name=strrep(filename,ee,['-cluster' get(handles.nclust,'String')]);
%    epsprint(1,name,1);
   exportpng1(1,[name '.png']);
end



% --- Executes on button press in checkall.
function checkall_Callback(hObject, eventdata, handles)
val=1.0;
set(handles.checkbox2,'Value',val);
set(handles.checkbox3,'Value',val);
set(handles.checkbox4,'Value',val);
set(handles.checkbox5,'Value',val);
set(handles.checkbox6,'Value',val);
set(handles.checkbox7,'Value',val);
set(handles.checkbox8,'Value',val);
set(handles.checkbox9,'Value',val);
set(handles.checkbox10,'Value',val);
set(handles.checkbox11,'Value',val);
set(handles.checkbox12,'Value',val);
set(handles.checkbox13,'Value',val);
set(handles.checkbox14,'Value',val);
set(handles.checkbox15,'Value',val);
set(handles.checkbox16,'Value',val);
set(handles.checkbox17,'Value',val);
set(handles.checkbox18,'Value',val);
set(handles.checkbox19,'Value',val);
set(handles.checkbox20,'Value',val);
set(handles.checkbox21,'Value',val);
set(handles.checkbox22,'Value',val);
set(handles.checkbox23,'Value',val);
set(handles.checkbox24,'Value',val);
set(handles.checkbox25,'Value',val);
set(handles.checkbox26,'Value',val);
set(handles.checkbox27,'Value',val);
set(handles.checkbox28,'Value',val);
set(handles.checkbox29,'Value',val);
set(handles.checkbox30,'Value',val);

% --- Executes on button press in checknone.
function checknone_Callback(hObject, eventdata, handles)
val=0.0;
set(handles.checkbox2,'Value',val);
set(handles.checkbox3,'Value',val);
set(handles.checkbox4,'Value',val);
set(handles.checkbox5,'Value',val);
set(handles.checkbox6,'Value',val);
set(handles.checkbox7,'Value',val);
set(handles.checkbox8,'Value',val);
set(handles.checkbox9,'Value',val);
set(handles.checkbox10,'Value',val);
set(handles.checkbox11,'Value',val);
set(handles.checkbox12,'Value',val);
set(handles.checkbox13,'Value',val);
set(handles.checkbox14,'Value',val);
set(handles.checkbox15,'Value',val);
set(handles.checkbox16,'Value',val);
set(handles.checkbox17,'Value',val);
set(handles.checkbox18,'Value',val);
set(handles.checkbox19,'Value',val);
set(handles.checkbox20,'Value',val);
set(handles.checkbox21,'Value',val);
set(handles.checkbox22,'Value',val);
set(handles.checkbox23,'Value',val);
set(handles.checkbox24,'Value',val);
set(handles.checkbox25,'Value',val);
set(handles.checkbox26,'Value',val);
set(handles.checkbox27,'Value',val);
set(handles.checkbox28,'Value',val);
set(handles.checkbox29,'Value',val);
set(handles.checkbox30,'Value',val);

% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
if ishandle(1), close(1); end
delete(gcbf);


% --- Executes during object creation, after setting all properties.
function distmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in distmethod.
function distmethod_Callback(hObject, eventdata, handles)
% hObject    handle to distmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns distmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from distmethod


% --- Executes during object creation, after setting all properties.
function linkmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to linkmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in linkmethod.
function linkmethod_Callback(hObject, eventdata, handles)
% hObject    handle to linkmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns linkmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from linkmethod


