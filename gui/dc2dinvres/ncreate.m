function varargout = ncreate(varargin)
% NCREATE Application M-file for ncreate.fig
%    FIG = NCREATE launch ncreate GUI.
%    NCREATE('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 30-Oct-2004 21:48:40

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    iconify(fig);
    global N
    if ~isempty(N),
        stat=[num2str(size(N.elec,1)) 'electrodes, ' num2str(length(N.a)) 'data'];
        set(handles.status,'String',stat);
        set(handles.nel,'String',num2str(size(N.elec,1)));
        set(handles.spacing,'String',num2str(min(diff(N.elec(:,1)))));
        set(handles.first,'String',num2str(min(N.elec(:,1))));
    end

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end

% --------------------------------------------------------------------
function varargout = arr_Callback(h, eventdata, handles, varargin)
arr=get(handles.arr,'Value');
onoff='Off';
if arr>9, onoff='On'; end
set(handles.bpostext,'Visible',onoff);
set(handles.bpos,'Visible',onoff);
set(handles.bdeptext,'Visible',onoff);
set(handles.bdeps,'Visible',onoff);
set(handles.bsptext,'Visible',onoff);
set(handles.bspacing,'Visible',onoff);
set(handles.tdip,'Visible',onoff);
set(handles.rdip,'Visible',onoff);
set(handles.text21,'Visible',onoff);
onoff='Off';
enl=get(handles.enlarge,'Value');
if ismember(arr,[3 6 7 8 9]), onoff='On'; end
set(handles.enlarge,'Visible',onoff),
if ~enl, onoff='Off'; end
set(handles.every,'Visible',onoff);
set(handles.shift,'Visible',onoff);
set(handles.text27,'Visible',onoff);
set(handles.text28,'Visible',onoff);

% --------------------------------------------------------------------
function varargout = nel_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = el1_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = del_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = nmin_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = nmax_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = bpos_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = bdeps_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = bspacing_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function add_Callback(hObject, eventdata, handles)
global N
nn=N;
ncreate('create_Callback',gcbo,[],guidata(gcbo));
N=combdata2d(nn,N);
stat=[num2str(size(N.elec,1)) 'elec, ' num2str(length(N.a)) 'data, err='...
        num2str(min(N.err)*100,'%.1f') '%-' num2str(max(N.err)*100,'%.1f') '%'];
set(handles.status,'String',stat);

% --------------------------------------------------------------------
function varargout = create_Callback(h, eventdata, handles, varargin)
global N
eel=str2num(get(handles.first,'String')); % first electrode
del=str2num(get(handles.spacing,'String')); % electrode spacing
nel=str2num(get(handles.nel,'String')); % number of electrodes
N.elec=(0:nel-1)'*del+eel;N.elec(1,2)=0;
arr=get(handles.arr,'Value'); % arrangement
nmin=str2num(get(handles.nmin,'String')); % minimum separation
nmax=str2num(get(handles.nmax,'String')); % maximum separation
N.a=[];N.b=[];N.m=[];N.n=[];
enl=get(handles.enlarge,'Value');
every=str2num(get(handles.every,'String'));
shift=str2num(get(handles.shift,'String'));
if arr<10, % surface measurements
    for n=nmin:nmax,
        df=1;ddf=1; % dipole enlargement for high separations
        if enl&ismember(arr,[3 6 7 8 9]), % DD,PD(f/r),CCPPC,SCHLUM
%             df=fix(n/(every+1))+1;
%             ddf=fix(n/(shift+1))+1;
            if every>0, df=ceil(n/every); end
            if shift==0, ddf=1; else ddf=ceil(n/shift); end
        end
        first=(1:ddf:nel)';
        abmn=zeros(length(first),4);
        abmn(:,1)=first;
        switch arr,
            case 1, %Wenner
                abmn(:,3)=abmn(:,1)+n;
                abmn(:,4)=abmn(:,1)+2*n;
                abmn(:,2)=abmn(:,1)+3*n;
            case 2, %Pole-Pole
                abmn(:,3)=abmn(:,1)+n*df;
            case {3,9}, %Dipole-Dipole and CC-PP-C
                abmn(:,2)=abmn(:,1)+df;
                abmn(:,3)=abmn(:,2)+n;
                abmn(:,4)=abmn(:,3)+df;
            case 4, %Wenner-beta
                abmn(:,2)=abmn(:,1)+n;
                abmn(:,3)=abmn(:,1)+2*n;
                abmn(:,4)=abmn(:,1)+3*n;
            case 5, %Wenner-gamma
                abmn(:,3)=abmn(:,1)+n;
                abmn(:,2)=abmn(:,1)+2*n;
                abmn(:,4)=abmn(:,1)+3*n;
            case 6, %Pole-Dipole
                abmn(:,3)=abmn(:,1)+n;
                abmn(:,4)=abmn(:,3)+df;
            case 8, %Pole-Dipole reverse
                abmn(:,1)=nel+1-first;
                abmn(:,3)=abmn(:,1)-n;
                abmn(:,4)=abmn(:,3)-df;
                abmn(find(abmn(:,4)<1),4)=-1;
            case 7, %Schlumberger
                abmn(:,3)=abmn(:,1)+n;
                abmn(:,4)=abmn(:,3)+df;
                abmn(:,2)=abmn(:,4)+n;
        end % switch
        fi=find(max(abmn')>nel);
        abmn(fi,:)=[];
        fi=find(min(abmn')<0);
        abmn(fi,:)=[];
        if ~isempty(abmn),
            N.a=[N.a;abmn(:,1)];
            N.b=[N.b;abmn(:,2)];
            N.m=[N.m;abmn(:,3)];
            N.n=[N.n;abmn(:,4)];
        end
    end % n loop
    if arr==9, % circulating dipole
        nradd=nel-3;
        N.a=[N.a;ones(nradd,1)];
        N.b=[N.b;ones(nradd,1)*nel];
        N.m=[N.m;(1:nradd)'+1];
        N.n=[N.n;(1:nradd)'+2];
    end
else
    %set(handles.status,'String','Boreholes not yet implemented!');
    bpos=str2num(get(handles.bpos,'String'));
    bdeps=str2num(get(handles.bdeps,'String'));
    bspac=str2num(get(handles.bspacing,'String'));
    tdip=get(handles.tdip,'Value');
    rdip=get(handles.rdip,'Value');
    N=createsxhdata(nel*(arr<11),eel,del,bpos,fix(bdeps/bspac),rdip,tdip);
    N=createsxhdata(nel*(arr<11),eel,bspac,bpos,fix(bdeps/bspac),rdip,tdip);
end
res=str2num(get(handles.reslevel,'String'));
N.r=ones(size(N.a))*res;
if get(handles.ip,'Value'), N.ip=zeros(size(N.a)); end
N.k=getkonf(N);
proz=str2num(get(handles.proz,'String'));
umin=str2num(get(handles.umin,'String'));
current=str2num(get(handles.current,'String'));
N.err=estimateerror(N,proz,umin/1000,current/1000);
elvar=str2num(get(handles.elvar,'String'));
if elvar>0,
    maxerr=electrodeerr2d(N,elvar);
    N.err=N.err+maxerr;
end
stat=[num2str(size(N.elec,1)) 'elec, ' num2str(length(N.a)) 'data, err='...
        num2str(min(N.err)*100,'%.1f') '%-' num2str(max(N.err)*100,'%.1f') '%'];
set(handles.status,'String',stat);

% --------------------------------------------------------------------
function varargout = cancel_Callback(h, eventdata, handles, varargin)
delete(gcbf);  %close figure

% --- Executes on button press in show.
function show_Callback(hObject, eventdata, handles)
global N
if length(N.a)==0,
    ncreate('create_Callback',gcbo,[],guidata(gcbo))
end
showdata2d(N,N.err*100,struct('log',1,'elec',1));
title('Error estimate in %');


% --- Executes on button press in tdip.
function tdip_Callback(hObject, eventdata, handles)

% --- Executes on button press in rdip.
function rdip_Callback(hObject, eventdata, handles)

% --- Executes on button press in enlarge.
function enlarge_Callback(hObject, eventdata, handles)
enl=get(handles.enlarge,'Value');
onoff='Off';
if enl, 
    onoff='On';
    set(handles.every,'String','4');
    set(handles.shift,'String','4');
end
set(handles.every,'Visible',onoff);
set(handles.shift,'Visible',onoff);
set(handles.text27,'Visible',onoff);
set(handles.text28,'Visible',onoff);

% --- Executes during object creation, after setting all properties.
function proz_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function proz_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function umin_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function umin_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function current_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function current_Callback(hObject, eventdata, handles)

function nmin_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function nmax_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function nel_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function el1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function del_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function arr_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function reslevel_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function reslevel_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function elvar_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function elvar_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function every_CreateFcn(hObject, eventdata, handles)
% hObject    handle to every (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function every_Callback(hObject, eventdata, handles)
% hObject    handle to every (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of every as text
%        str2double(get(hObject,'String')) returns contents of every as a double


% --- Executes during object creation, after setting all properties.
function shift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function shift_Callback(hObject, eventdata, handles)


% --- Executes on button press in ip.
function ip_Callback(hObject, eventdata, handles)

