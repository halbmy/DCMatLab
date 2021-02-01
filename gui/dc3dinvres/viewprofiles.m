function varargout = viewprofiles(varargin)

% Last Modified by GUIDE v2.5 12-Aug-2006 00:43:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @viewprofiles_OpeningFcn, ...
                   'gui_OutputFcn',  @viewprofiles_OutputFcn, ...
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

% --- Executes just before viewprofiles is made visible.
function viewprofiles_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for viewprofiles
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using viewprofiles.
global N field
set(handles.profile,'String',N.names);
if ~isequal(size(field),size(N.r)), field=N.r; end
viewprofiles('profile_Callback',hObject,[],handles); 

% UIWAIT makes viewprofiles wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = viewprofiles_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
delete(handles.figure1)

% --- Executes during object creation, after setting all properties.
function field_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in popupmenu3.
function field_Callback(hObject, eventdata, handles)
global N Model field
switch get(hObject,'Value'),
    case 2, field=N.r;
    case 3, field=Model.R;
    case 4, field=(1-Model.R./N.r)*100;
    case 5, field=N.err;
    case 6, field=N.i;
    case 7, field=N.u;
    case 8, field=N.ip;
end
viewprofiles('profile_Callback',gcbo,[],guidata(gcbo)); 


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
prof=get(handles.profile,'Value');
if prof<length(get(handles.profile,'String')),
    set(handles.profile,'Value',prof+1);
end
viewprofiles('profile_Callback',gcbo,[],guidata(gcbo)); 

% --- Executes on button press in previous.
function previous_Callback(hObject, eventdata, handles)
prof=get(handles.profile,'Value');
if prof>1,
    set(handles.profile,'Value',prof-1);
end
viewprofiles('profile_Callback',gcbo,[],guidata(gcbo)); 

% --- Executes during object creation, after setting all properties.
function profile_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in profile.
function profile_Callback(hObject, eventdata, handles)
global N field MAL
prof=get(handles.profile,'Value');
toshow=field(N.nr{prof});
mal.clog=(min(toshow)>0);
mal.cauto=1;
if get(handles.fix,'Value'),
    mal.cauto=0;
    if min(field)<0,
        mal.clog=0;mm=max(abs(field));mal.cmin=-mm;mal.cmax=mm;
    else
        mal.cmin=min(field);mal.cmax=max(field);
    end
end
set(handles.status,'String',sprintf('%d data with %d electrodes',length(N.zweid{prof}.r),size(N.zweid{prof}.elec,1)));
axes(handles.axes1);
showdata2d(N.zweid{prof},toshow,mal);
im=get(handles.axes1,'Children');
set(im,'ButtonDownFcn','viewprofiles(''malfeld_ButtonDownFcn'',gcbo,[],guidata(gcbo))') 
axes(handles.axes2);
xl='x/m';yl='y/m';
px=1;py=2;xy='XY';
if isfield(MAL,'xy')&&(MAL.xy>0),
    du=xl;xl=yl;yl=du;xy='YX';
    du=px;px=py;py=du;
end
for i=1:length(N.points),
    po=N.points{i};
    col='b';if i==prof, col='r'; end
    plot(po(px:2:end),po(py:2:end),[col '-']);
    if i==1, hold on; end
    plot(po(px(1)),po(py(1)),[col '.']);
end
hold off

axis equal tight
if isfield(MAL,'xdir')&&(MAL.xdir==1), set(gca,[xy(1) 'Dir'],'reverse'); end
if isfield(MAL,'ydir')&&(MAL.ydir==1), set(gca,[xy(2) 'Dir'],'reverse'); end
xtl=get(gca,'XTickLabel');if ischar(xtl), xtl=cellstr(xtl); end
xtl{end-1}=xl;set(gca,'XTickLabel',xtl);
ytl=get(gca,'YTickLabel');if ischar(ytl), ytl=cellstr(ytl); end
ytl{end-1}=yl;set(gca,'YTickLabel',ytl);

function varargout = malfeld_ButtonDownFcn(h, eventdata, handles, varargin) 
global N field
cp=get(handles.axes1,'CurrentPoint');
prof=get(handles.profile,'Value');
NN=N.zweid{prof};
% NN.r=field(N.nr{prof});
[mids,konfs,ii,kk]=midkonf2d(NN); 
nk=round(cp(1,2));
[aa,nx]=min(abs(mids-cp(1,1)));
nr=find((kk==nk)&(ii==nx));
if ~isempty(nr),
    if length(nr)>1, nr=nr(end); end
    mess=sprintf('Nr %d: Model.R=%.1f',nr,NN.r(nr));
    if isfield(NN,'ip')&&(length(NN.ip)>=nr), 
        mess=sprintf('%s IP=%.2fmrad',mess,NN.ip(nr)); end  
    if isfield(NN,'i')&&(length(NN.i)>=nr), 
        mess=sprintf('%s I=%.1fmA',mess,NN.i(nr)*1000); end  
    if isfield(NN,'u')&&(length(NN.u)>=nr), 
        mess=sprintf('%s U=%.1fmV',mess,NN.u(nr)*1000); end  
    if isfield(NN,'err')&&(length(NN.err)>=nr), 
        mess=sprintf('%s Error=%.1f%%',mess,NN.err(nr)*100); end  
    mess=sprintf('%s\na=%d (%.1fm)',mess,NN.a(nr),NN.elec(NN.a(nr),1));
    if NN.b(nr)>0, mess=sprintf('%s  b=%d (%.1fm)',mess,NN.b(nr),NN.elec(NN.b(nr),1)); end
    mess=sprintf('%s  m=%d (%.1fm)',mess,NN.m(nr),NN.elec(NN.m(nr),1));
    if NN.n(nr)>0, mess=sprintf('%s  n=%d (%.1fm)',mess,NN.n(nr),NN.elec(NN.n(nr),1)); end
    if strcmp(questdlg(mess,'Delete Datum?','Yes','No','No'),'Yes'),
        NN.a(nr)=[];NN.b(nr)=[];NN.m(nr)=[];NN.n(nr)=[];
        global Model INV S RMS
        NN.r(nr)=[];NN.k(nr)=[];
        if isfield(NN,'err')&&(length(NN.err)>=nr), NN.err(nr)=[]; end
        if isfield(NN,'ip')&&(length(NN.ip)>=nr), NN.ip(nr)=[]; end
        if isfield(NN,'i')&&(length(NN.i)>=nr), NN.i(nr)=[]; end
        if isfield(NN,'u')&&(length(NN.u)>=nr), NN.u(nr)=[]; end
        N.zweid{prof}=NN;
        nnr=N.nr{prof}(nr);
        map=[0:nnr-1 0 nnr:length(N.r)];
        for i=1:length(N.nr), N.nr{i}=map(N.nr{i}+1); end
        N.nr{prof}(nr)=[];
        N.a(nnr)=[];N.b(nnr)=[];N.m(nnr)=[];N.n(nnr)=[];
        N.r(nnr)=[];N.k(nnr)=[];Model.R(nnr)=[];field(nnr)=[];
        if isfield(N,'err')&&(length(N.err)>=nnr), N.err(nnr)=[]; end
        if isfield(N,'ip')&&(length(N.ip)>=nnr), N.ip(nnr)=[]; end
        if isfield(N,'i')&&(length(N.i)>=nnr), N.i(nnr)=[]; end
        if isfield(N,'u')&&(length(N.u)>=nnr), N.u(nnr)=[]; end
        if size(S,1)>=nnr, S(nnr,:)=[]; end
%          message(sprintf('Deleting Datum Number %d from data set. new Chi^2=%.1f RMS=%.1f',nnr,CHIQ(end),RMS(end))); 
%         axes(handles.axes1);showdata2d(NN);      
        viewprofiles('profile_Callback',gcbo,[],guidata(gcbo)); 
    end
end

function varargout = viewprofiles_KeyPressFcn(h, eventdata, handles, varargin)
aa=get(gcf,'CurrentCharacter');
switch upper(aa), 
    case 'N',
        viewprofiles('next_Callback',gcbo,[],guidata(gcbo)); 
    case 'P',
        viewprofiles('previous_Callback',gcbo,[],guidata(gcbo)); 
    case 'X',
        viewprofiles('CloseMenuItem_Callback',gcbo,[],guidata(gcbo)); 
end

% --- Executes on button press in fix.
function fix_Callback(hObject, eventdata, handles)
viewprofiles('profile_Callback',gcbo,[],guidata(gcbo)); 


% --- Executes when figure1 window is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


