function varargout = s_inv(varargin)
% S_INV Application M-file for s_inv.fig
%    FIG = S_INV launch s_inv GUI.
%    S_INV('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 11-Aug-2006 10:24:28

if nargin == 0  % LAUNCH GUI
	fig = openfig(mfilename,'reuse');
	% Use system color scheme for figure:  
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    iconify(fig);
    global INV
    if isfield(INV,'mico'), set(handles.mico,'String',num2str(INV.mico)); end
    if isfield(INV,'lam'), 
      set(handles.regutext,'String',num2str(INV.lam,4)); 
      rr=-0.1;
      if INV.lam>=1, rr=log10(INV.lam); end
      ma=get(handles.regu,'Max');if rr>ma, rr=ma; end
      set(handles.regu,'Value',rr);
    end
    if isfield(INV,'lolo'), 
        set(handles.logrhoa,'Value',(INV.lolo~=0));
        set(handles.rhoa,'Value',(INV.lolo==0)); 
        if INV.lolo==0, 
            set(handles.lbound,'Visible','Off');
            set(handles.ubound,'Visible','Off');
            set(handles.text23,'Visible','Off');
            set(handles.text17,'Visible','Off');            
            set(handles.lboundtext,'Visible','Off'); end
    end
    if isfield(INV,'redu'), set(handles.redu,'Value',INV.redu+1); end
    if isfield(INV,'weight'), set(handles.weight,'Value',INV.weight+1); end
    if isfield(INV,'mitschicht'), set(handles.mitschicht,'Value',INV.mitschicht); end
    if isfield(INV,'const'), set(handles.const,'Value',INV.const); end
    if isfield(INV,'lbound'), set(handles.lbound,'String',num2str(INV.lbound)); end
    if isfield(INV,'ubound'), set(handles.ubound,'String',num2str(INV.ubound)); end
    if isfield(INV,'rbzfak'), set(handles.zfak,'String',num2str(INV.rbzfak)); end
    if isfield(INV,'linesearch'), set(handles.linesearch,'Value',INV.linesearch>0); end
    if isfield(INV,'glob'), set(handles.glob,'Value',INV.glob>0); end
    if isfield(INV,'robust'), set(handles.robustdata,'Value',INV.robust>0); end
    if isfield(INV,'blocky'), set(handles.blockymodel,'Value',INV.blocky>0); end
    if isfield(INV,'spsens')&&(INV.spsens>0),
        set(handles.issparse,'Value',1);
        set(handles.text19,'Enable','On');
        set(handles.spsens,'Enable','On','String',num2str(INV.spsens));
    end
    if isfield(INV,'auto'),
        set(handles.auto,'Value',INV.auto+1);
        if INV.auto==2,
            set(handles.text16,'String','lambda'); 
            set(handles.const,'Visible','Off');
        else
            set(handles.text16,'String','min. lambda'); 
            set(handles.const,'Visible','On');
        end
    end
    if isfield(INV,'method'), 
        set(handles.method,'Value',INV.method+1);
        s_inv('method_Callback',fig,[],guidata(fig));
    end
	if nargout > 0
		varargout{1} = fig;
	end
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end
end

% --------------------------------------------------------------------
function varargout = rhoa_Callback(h, eventdata, handles, varargin)
set(handles.logrhoa,'Value',0);
set(handles.lbound,'Visible','Off');
set(handles.ubound,'Visible','Off');
set(handles.text17,'Visible','Off');
set(handles.text23,'Visible','Off');
set(handles.lboundtext,'Visible','Off');

% --------------------------------------------------------------------
function varargout = logrhoa_Callback(h, eventdata, handles, varargin)
set(handles.rhoa,'Value',0);
set(handles.lbound,'Visible','On');
set(handles.ubound,'Visible','On');
set(handles.text17,'Visible','On');
set(handles.text23,'Visible','On');
set(handles.lboundtext,'Visible','On');

function varargout = redu_Callback(h, eventdata, handles, varargin)

function varargout = mitschicht_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = inv_ok_Callback(h, eventdata, handles, varargin)
global INV
INV.mico=str2num(get(handles.mico,'String'));
INV.method=get(handles.method,'Value')-1;
INV.mitschicht=get(handles.mitschicht,'Value');
INV.redu=get(handles.redu,'Value')-1;
INV.weight=get(handles.weight,'Value')-1;
INV.lolo=get(handles.logrhoa,'Value');
INV.lam=str2num(get(handles.regutext,'String'));
INV.auto=get(handles.auto,'Value')-1;
INV.const=get(handles.const,'Value');
INV.lbound=str2num(get(handles.lbound,'String'));
INV.ubound=str2num(get(handles.ubound,'String'));
INV.rbzfak=str2num(get(handles.zfak,'String'));
INV.linesearch=get(handles.linesearch,'Value');
INV.glob=get(handles.glob,'Value');
INV.robust=get(handles.robustdata,'Value');
INV.blocky=get(handles.blockymodel,'Value');
INV.spsens=str2num(get(handles.spsens,'String'))*get(handles.issparse,'Value');
%if INV.method==0, INV.lam=0; end
delete(gcbf);

% --------------------------------------------------------------------
function varargout = inv_cancel_Callback(h, eventdata, handles, varargin)
delete(gcbf);

function varargout = regu_Callback(h, eventdata, handles, varargin)
rr=get(handles.regu,'Value');
lam=0;
if rr>=0, lam=10^(rr); end
set(handles.regutext,'String',num2str(lam,3));

function varargout = regutext_Callback(h, eventdata, handles, varargin)
mi=get(handles.regu,'Min');
ma=get(handles.regu,'Max');
rr=log10(str2num(get(handles.regutext,'String')));
if (rr>=mi)&(rr<=ma), set(handles.regu,'Value',rr); end
% if isnumeric(val) & length(val)==1 & val >= get(handles.regu,'Min') & ...
%     val <= get(handles.regu,'Max')
%     set(handles.regu,'Value',val);
% else
%     set(handles.regutext,'String',num2str(get(handles.regu),'Value'));
% end

function varargout = method_Callback(h, eventdata, handles, varargin)
meth=get(handles.method,'Value');
if meth>1, status='Off'; else status='On'; end
status1=status;if ismember(meth,[2 4]), status1='On'; end
set(handles.weight,'Visible',status);
set(handles.redu,'Visible',status);
set(handles.mitschicht,'Visible',status);
set(handles.const,'Visible',status1);
set(handles.lbound,'Visible',status1);
set(handles.ubound,'Visible',status1);
set(handles.text17,'Visible',status1);
set(handles.text23,'Visible',status1);
set(handles.lboundtext,'Visible',status1);
set(handles.text16,'Visible',status);
set(handles.text13,'Visible',status);
set(handles.text2,'Visible',status);
set(handles.text12,'Visible',status);
set(handles.auto,'Visible',status1);
set(handles.mico,'Visible',status);
set(handles.regu,'Visible',status1);
set(handles.regutext,'Visible',status1);

% --------------------------------------------------------------------
function varargout = auto_Callback(h, eventdata, handles, varargin)
if get(handles.auto,'Value')==3, %manual choosing
    set(handles.text16,'String','lambda'); 
    set(handles.const,'Visible','Off');
    set(handles.regutext,'String','30');
    s_inv('regutext_Callback',gcbo,[],guidata(gcbo));
else
    set(handles.text16,'String','min. lambda'); 
    set(handles.const,'Visible','On');
    set(handles.regutext,'String','1');
    s_inv('regutext_Callback',gcbo,[],guidata(gcbo));
end

% --- Executes during object creation, after setting all properties.
function weight_CreateFcn(hObject, eventdata, handles)

% --- Executes on selection change in weight.
function weight_Callback(hObject, eventdata, handles)

% --- Executes on button press in const.
function const_Callback(hObject, eventdata, handles)


% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
ibrowse(['file://' pwd filesep 'doc' filesep 'options.html#inv'])

% --- Executes during object creation, after setting all properties.
function lbound_CreateFcn(hObject, eventdata, handles)
% if ispc
%     set(hObject,'BackgroundColor','white');
% else
%     set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
% end

function lbound_Callback(hObject, eventdata, handles)


% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
aa=get(handles.figure1,'CurrentCharacter');
if double(aa)==27, s_inv('inv_cancel_Callback',gcbo,[],guidata(gcbo)); end
if double(aa)==13, s_inv('inv_ok_Callback',gcbo,[],guidata(gcbo)); end


% --- Executes on button press in linesearch.
function linesearch_Callback(hObject, eventdata, handles)


% --- Executes on button press in glob.
function glob_Callback(hObject, eventdata, handles)


% --- Executes on button press in issparse.
function issparse_Callback(hObject, eventdata, handles)
global INV
onoff='Off';
if get(handles.issparse,'Value')>0,
    onoff='On';
    if isfield(INV,'spsens')&&(INV.spsens>0),
        set(handles.spsens,'String',num2str(INV.spsens));
    end
end
set(handles.text19,'Enable',onoff);
set(handles.spsens,'Enable',onoff);

% --- Executes during object creation, after setting all properties.
function spsens_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spsens (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function spsens_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function zfak_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zfak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function zfak_Callback(hObject, eventdata, handles)
% hObject    handle to zfak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zfak as text
%        str2double(get(hObject,'String')) returns contents of zfak as a double


% --- Executes during object creation, after setting all properties.
function ubound_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ubound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function ubound_Callback(hObject, eventdata, handles)
% hObject    handle to ubound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in robustdata.
function robustdata_Callback(hObject, eventdata, handles)


% --- Executes on button press in blockymodel.
function blockymodel_Callback(hObject, eventdata, handles)
