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
if nargin & isstr(varargin{1})
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
    ss{end+1}='Reciprocity';
    set(handles.field,'String',ss);
end
if isfield(N,'u')&&(length(N.u)==length(N.a)),
    ss=get(handles.field,'String');
    ss{end+1}='Voltage';
    set(handles.field,'String',ss);
end
if isfield(N,'i')&&(length(N.i)==length(N.a)),
    ss=get(handles.field,'String');
    ss{end+1}='Current';
    set(handles.field,'String',ss);
end
ss=sprintf('%d data points present',length(N.r));
set(handles.status,'String',ss);
if ~isfield(N,'zweid')||(length(N.zweid)<15),
    ndelete('show_Callback',hObject,[],guidata(hObject));
end
ndelete('status_Callback',hObject,[],guidata(hObject));

% UIWAIT makes ndelete wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = ndelete_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
%varargout{1} = handles.output;

% --- Executes on button press in delete.
function delete_Callback(hObject, eventdata, handles)
global N fdel Model S
if ~isempty(fdel),
    ndata=length(N.r);
    map=zeros(ndata,1);
    map(setxor(1:ndata,fdel))=1:ndata-length(fdel);
    N.a(fdel)=[];N.b(fdel)=[];N.m(fdel)=[];N.n(fdel)=[];
    N.r(fdel)=[];N.k(fdel)=[];
    if isfield(N,'rez')&&(length(N.rez)>=max(fdel)), N.rez(fdel)=[]; end
    if isfield(N,'err')&&(length(N.err)>=max(fdel)), N.err(fdel)=[]; end
    if isfield(N,'ip')&&(length(N.ip)>=max(fdel)), N.ip(fdel)=[]; end
    if isfield(N,'u')&&(length(N.u)>=max(fdel)), N.u(fdel)=[]; end
    if isfield(N,'i')&&(length(N.i)>=max(fdel)), N.i(fdel)=[]; end
    if isfield(N,'rho')&&(length(N.rho)>=max(fdel)), N.rho(fdel)=[]; end
    if length(Model.R)>=max(fdel), Model.R(fdel)=[]; end
    if size(S,1)>=max(fdel), S(fdel,:)=[]; end
%     map=[0;map]; % wherever that comes from... it's bad!
    if isfield(N,'zweid'),
        for n=1:length(N.zweid),
            nn=N.zweid{n};
            N.nr{n}=map(N.nr{n});
            fi=find(N.nr{n}==0);
            if ~isempty(fi),
                N.nr{n}(fi)=[];
                %             nn.a=map(nn.a+1)-1;nn.b=map(nn.b+1)-1;
                %             nn.m=map(nn.m+1)-1;nn.n=map(nn.n+1)-1;
                nn.a(fi)=[];nn.b(fi)=[];
                nn.m(fi)=[];nn.n(fi)=[];
                nn.r(fi)=[];nn.k(fi)=[];
                if isfield(nn,'err')&&(length(nn.err)>=max(fi)), nn.err=[]; end
                if isfield(nn,'ip')&&(length(nn.ip)>=max(fi)), nn.ip=[]; end
                N.zweid{n}=nn;
            end
        end
    end
    
    if ~isfield(N,'zweid')||(length(N.zweid)<15),
        ndelete('show_Callback',gcbo,[],guidata(gcbo));
    end
    ss=sprintf('Deleted %d data',length(fdel));
    set(handles.status,'String',ss);
else
    set(handles.status,'String','Nothing to delete!');
end

% --- Executes on button press in show.
function show_Callback(hObject, eventdata, handles)
global N MAL Model
fie=get(handles.field,'Value');
mal=MAL;mal.cauto=1;mal.clog=1;
feld=N.r;figtitle='Data';
if fie==2, % error
    feld=N.err*100;figtitle='Error'; 
    mal.clog=0;
end
if (fie==6)&&isfield(N,'ip')&&(length(N.ip)==length(N.r)),
    feld=N.ip;figtitle='IP data';
    mal.clog=0;
end
if fie==7,
    feld=(1-Model.R./N.r)*100;figtitle='Misfit';
    mal.clog=0;mal.cmap=2;
end
if fie==8,
    feld=N.k;figtitle='Geom. factor';
    mal.clog=0;mal.cmap=2;
end
if fie>=9,
    ss=get(handles.field,'String');
    figtitle=ss{fie};
    if strcmp(figtitle,'Reciprocity'),
        feld=N.rez*100;
        mal.clog=0;mal.cmap=2;
    end
    if strcmp(figtitle,'Voltage'),
        feld=N.u;
        mal.clog=0;
    end
    if strcmp(figtitle,'Current'),
        feld=abs(N.i);
        mal.clog=1;
    end
end
set(figure(2),'MenuBar','none','NumberTitle','off','Name','Data');
iconify(2);
if isfield(N,'zweid')||isfield(N,'eind'), plotprofiles(N,feld,mal); 
else hist(feld,30); end

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
if fie==1, ss='Om'; end
if ismember(fie,[2 7 8]), ss='%'; end
if ismember(fie,3:5), ss='m'; end
set(handles.unit,'String',ss);
ndelete('status_Callback',gcbo,[],guidata(gcbo));
if ~isfield(N,'zweid')||(length(N.zweid)<15),
    ndelete('show_Callback',gcbo,[],guidata(gcbo));
end

function status_Callback(hObject, eventdata, handles)
global N MAL fdel zeig Model
rel=get(handles.relation,'Value');
fie=get(handles.field,'Value');
val=str2num(get(handles.wert,'String'));
feld=N.r;lolo=min(N.r>0); % resistivity
if fie==2, feld=N.err*100; end % error
if fie==7, feld=(Model.R./N.r-1)*100; end % misfit
if fie==8, feld=N.k; end
if fie>=9,
    ss=get(handles.field,'String');
    figtitle=ss{fie};
    if strcmp(figtitle,'Reciprocity'), feld=abs(N.rez)*100; end
    if strcmp(figtitle,'Voltage'), feld=N.u; end
    if strcmp(figtitle,'Current'), feld=abs(N.i); end
end
if (fie==6)&&isfield(N,'ip')&&(length(N.ip)==length(N.r)), feld=N.ip; end
if ismember(fie,3:5), % x/y/z position
    abmn=zeros(length(N.r),4);
    abmn(:,1)=N.elec(N.a,fie-2);
    abmn(:,2)=abmn(:,1);
    fi=find(N.b);
    abmn(fi,2)=N.elec(N.b(fi),fie-2);
    abmn(:,3)=N.elec(N.m,fie-2);
    abmn(:,4)=abmn(:,3);
    fi=find(N.n);
    abmn(fi,4)=N.elec(N.n(fi),fie-2);
    if rel==1, feld=max(abmn,[],2); end
    if rel==2, feld=min(abmn,[],2); end
    if rel==3, 
        feld=prod(abmn-val,2);
        val=0;
    end
end
if rel==1, % greater 
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
% ndelete('marked_Callback',hObject,[],handles);
mal=MAL;
mal.cauto=0;mal.cmin=0;mal.cmax=1;mal.clog=0;mal.clog=0;mal.cbar=0;
figure(2);
plotprofiles(N,zeig,mal);


% --- Executes on button press in stat.
function stat_Callback(hObject, eventdata, handles)
global N Model
fie=get(handles.field,'Value');
feld=N.r;lolo=min(N.r>0); % resistivity
if fie==2, feld=N.err*100;lolo=0; end % error
if (fie==6)&&isfield(N,'ip')&&(length(N.ip)==length(N.r)), feld=N.ip;lolo=0; end
if fie==7, feld=(Model.R./N.r-1)*100;lolo=0; end % misfit
if fie==8, feld=N.k;lolo=0; end
if fie>=9,
    ss=get(handles.field,'String');
    figtitle=ss{fie};lolo=0;
    if strcmp(figtitle,'Reciprocity'), feld=abs(N.rez)*100;lolo=1; end
    if strcmp(figtitle,'Voltage'), feld=N.u; end
    if strcmp(figtitle,'Current'), feld=abs(N.i);lolo=1; end
end
set(figure(2),'MenuBar','none','NumberTitle','off','Name','Data');
iconify(2);clf;
if lolo,
    hist(log10(feld),30);
    xl=rndig(10.^get(gca,'XTick'));
%     fi=find(xl>1);
%     xl(fi)=round(xl(fi)*10)/10;
%     fi=find(xl>20);
%     xl(fi)=round(xl(fi));
    set(gca,'XTicklabel',num2cell(xl));
else
    hist(feld,30);
end

% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
clear('global','fdel','zeig');
%uiresume(handles.figure1);
delete(gcbf);

