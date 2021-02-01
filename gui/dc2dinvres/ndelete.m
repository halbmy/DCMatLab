function varargout = ndelete(varargin)

% Last Modified by GUIDE v2.5 23-Jul-2003 15:44:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ndelete_OpeningFcn, ...
                   'gui_OutputFcn',  @ndelete_OutputFcn, ...
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

% --- Executes just before ndelete is made visible.
function ndelete_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for ndelete
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
iconify(hObject);
global N
if isfield(N,'rez')&&(length(N.rez)==length(N.a)),
    ss=get(handles.field,'String');
    ss{9}='Reciprocity';
    set(handles.field,'String',ss);
end
ss=sprintf('%d data points present',length(N.r));
set(handles.status,'String',ss);
set(handles.relation,'Value',2);
ndelete('show_Callback',hObject,[],guidata(hObject));

% UIWAIT makes ndelete wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = ndelete_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
%varargout{1} = handles.output;

% --- Executes on button press in delete.
function delete_Callback(hObject, eventdata, handles)
global N fdel Mod S eldel
if ~isempty(fdel),
%     if isfield(N,'err')&&(length(N.err)==length(N.a)), N.err(fdel)=[]; end
%     if isfield(N,'ip')&&(length(N.ip)==length(N.a)), N.ip(fdel)=[]; end
%     if isfield(N,'i')&&(length(N.i)==length(N.a)), N.i(fdel)=[]; end
%     if isfield(N,'u')&&(length(N.u)==length(N.a)), N.u(fdel)=[]; end
%     N.a(fdel)=[];N.b(fdel)=[];N.m(fdel)=[];N.n(fdel)=[];
%     N.r(fdel)=[];N.k(fdel)=[];
    fn=fieldnames(N);si=size(N.a);
    for i=1:length(fn), field=getfield(N,fn{i});
        if isequal(size(field),si), field(fdel)=[];N=setfield(N,fn{i},field); end; end
    if length(Mod.R)>=max(fdel), Mod.R(fdel)=[]; end
    if size(S,1)>=max(fdel), S(fdel,:)=[]; end
    ss=sprintf('Deleted %d data',length(fdel));
    if ~isempty(eldel),
        neldel=setxor(1:size(N.elec,1),eldel);
        map=ones(size(N.elec,1),1);    
        map(neldel)=(1:length(neldel))';
        N.a=map(N.a);
        N.m=map(N.m);
        fi=find(N.b);N.b(fi)=map(N.b(fi));
        fi=find(N.n);N.n(fi)=map(N.n(fi));
        N.elec(eldel,:)=[];
        ss=sprintf('%s , %d electrodes',ss,length(eldel));
    end
    fdel=[];
    set(handles.status,'String',ss);
    ndelete('show_Callback',gcbo,[],guidata(gcbo));
else
    set(handles.status,'String','Nothing to delete!');
end

% --- Executes on button press in show.
function show_Callback(hObject, eventdata, handles)
global N MAL
fie=get(handles.field,'Value');
mal=MAL;mal.cauto=1;
feld=N.r;
if fie==2, %err
    feld=N.err*100; 
    mal.clog=0;
end
if fie==4, %k-factor
    feld=abs(N.k); 
    mal.clog=1;
end
if (fie==5)&&isfield(N,'ip'), % IP
    feld=N.ip;mal.clog=0;
end
if (fie==6)&&isfield(N,'i'), % I
    feld=abs(N.i);mal.clog=1;
end
if (fie==7)&&isfield(N,'u'), % U
    feld=abs(N.u);mal.clog=1;
end
if (fie==8), % misfit
    global Mod
    feld=(1-Mod.R./N.r)*100;mal.cmap=2;
end
if (fie==9)&&isfield(N,'rez'), % Reciprocity
    feld=abs(N.rez)*100;mal.clog=0;
end
if isfield(mal,'clog'), mal.log=mal.clog; end
%figure(1);
showdata2d(N,feld,mal);

% --- Executes during object creation, after setting all properties.
function wert_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function wert_Callback(hObject, eventdata, handles)
ndelete('status_Callback',gcbo,[],guidata(gcbo));

% --- Executes during object creation, after setting all properties.
function relation_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in relation.
function relation_Callback(hObject, eventdata, handles)
ndelete('status_Callback',gcbo,[],guidata(gcbo));

% --- Executes during object creation, after setting all properties.
function field_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in field.
function field_Callback(hObject, eventdata, handles)
global N
fie=get(handles.field,'Value');
ss='';
if fie==1, ss='Om'; end %rhoa
if ismember(fie,[2 8 9]), ss='%'; end %err,misfit
if ismember(fie,[3 4]), ss='m'; end %pos,kfak
if fie==5, 
    if isfield(N,'ip'),
        ss='mrad'; 
    else
        set(handles.field,'Value',1);
        ndelete('field_Callback',gcbo,[],guidata(gcbo));
    end
end %IP
if fie==6, 
    if isfield(N,'i'),
        ss='A'; 
    else
        set(handles.field,'Value',1);
        ndelete('field_Callback',gcbo,[],guidata(gcbo));
    end
end %I
if fie==7, 
    if isfield(N,'u'),
        ss='V'; 
    else
        set(handles.field,'Value',1);
        ndelete('field_Callback',gcbo,[],guidata(gcbo));
    end
end %U
set(handles.unit,'String',ss);
ndelete('status_Callback',gcbo,[],guidata(gcbo));
ndelete('show_Callback',gcbo,[],guidata(gcbo));

function status_Callback(hObject, eventdata, handles)
global N MAL fdel zeig eldel Mod
rel=get(handles.relation,'Value');
fie=get(handles.field,'Value');
val=str2double(get(handles.wert,'String'));
feld=N.r; % resistivity
eldel=[];
if fie==2, feld=N.err*100; end % error
if fie==4, feld=abs(N.k); end % kfak
if fie==5, feld=N.ip; end
if fie==6, feld=abs(N.i); end
if fie==7, feld=abs(N.u); end
if fie==8, feld=(1-Mod.R./N.r)*100; end
if fie==9, feld=abs(N.rez)*100; end
if fie==3, % position
    abmn=zeros(length(N.r),4);
    abmn(:,1)=N.elec(N.a,1);
    abmn(:,2)=abmn(:,1);
    fi=find(N.b);
    abmn(fi,2)=N.elec(N.b(fi),1);
    abmn(:,3)=N.elec(N.m,1);
    abmn(:,4)=abmn(:,3);
    fi=find(N.n);
    abmn(fi,4)=N.elec(N.n(fi),1);
    if rel==1, % greater than
        feld=max(abmn,[],2); 
        eldel=find(N.elec(:,1)>val);
    end
    if rel==2, % less than
        feld=min(abmn,[],2); 
        eldel=find(N.elec(:,1)<val);
    end
    if rel==3, % equal
        feld=prod(abmn-val,2);
        eldel=find(N.elec(:,1)==val);
        val=0;
    end
end
if rel==1, % greater than 
    zeig=(feld>val);
end
if rel==2, % less than
    zeig=(feld<val);
end
if rel==3, % equal
    zeig=(feld==val);
end
fdel=find(zeig);
if ~isempty(fdel),
    ss=sprintf('Found %d data',sum(zeig));
else
    ss='No data found';
end
set(handles.status,'String',ss);

% --- Executes on button press in marked.
function marked_Callback(hObject, eventdata, handles)
global MAL N zeig
mal=MAL;
mal.cauto=0;mal.cmin=0;mal.cmax=1;mal.clog=0;mal.log=0;mal.cbar=0;
% figure(1);
showdata2d(N,zeig,mal);

% --- Executes on button press in stat.
function stat_Callback(hObject, eventdata, handles)
global N Mod
fie=get(handles.field,'Value');
feld=N.r;lolo=1; % resistivity
if fie==2, feld=N.err*100;lolo=0; end % error
if fie==4, feld=abs(N.k); end % k
if fie==5, feld=N.ip;lolo=0; end % 5
if fie==6, feld=abs(N.i); end % i
if fie==7, feld=abs(N.u); end % u
if fie==8, feld=(1-Mod.R./N.r)*100;lolo=0; end %misfit
if fie==9, feld=abs(N.rez)*100;lolo=0; end
% figure(1);
nn=max(30,fix(length(N.r)/30));
if lolo,
    hist(log10(feld),nn);
    xl=10.^get(gca,'XTick');
    fi=find(xl>1);
    xl(fi)=round(xl(fi)*10)/10;
    fi=find(xl>20);
    xl(fi)=round(xl(fi));
    set(gca,'XTicklabel',num2cell(xl));
else
    hist(feld,nn);
end

% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
global fdel zeig eldel
clear('global','fdel','zeig','eldel');
%uiresume(handles.figure1);
delete(gcbf);
