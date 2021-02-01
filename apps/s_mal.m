function varargout = s_mal(varargin)
% S_MAL Application M-file for s_mal.fig
%    FIG = S_MAL launch s_mal GUI.
%    S_MAL('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 19-Oct-2006 14:08:58

global MAL Mod
if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

% Use system color scheme for figure:
set(fig,'Color',get(0,'DefaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    iconify(fig);
   
    if isfield(MAL,'alpha'), set(handles.alpha,'Value',MAL.alpha); end
    if isfield(MAL,'xdir'), set(handles.xdir,'Value',MAL.xdir); end
    if isfield(MAL,'clog'), set(handles.log,'Value',MAL.clog); else set(handles.log,'Value',1); end
    if isfield(MAL,'cauto'), set(handles.cauto,'Value',MAL.cauto); else set(handles.cauto,'Value',1); end
    if isfield(MAL,'cmap'), set(handles.cmap,'Value',MAL.cmap+1); end
    if isfield(MAL,'cmin'), set(handles.cmin,'String',num2str(MAL.cmin)); end
    if isfield(MAL,'cmax'), set(handles.cmax,'String',num2str(MAL.cmax)); end
    if isfield(MAL,'creverse'), set(handles.creverse,'Value',MAL.creverse>0); end
    if isfield(MAL,'style'), set(handles.style,'Value',MAL.style+1); end
    if isfield(MAL,'showgrid'), set(handles.showgrid,'Value',MAL.showgrid+1); end
    if isfield(MAL,'high'), set(handles.high,'Value',MAL.high+1); end
    if (~isfield(MAL,'cauto'))||(MAL.cauto==1),
        set(handles.cmin,'Enable','off');
        set(handles.cmax,'Enable','off');
    end
    if isfield(MAL,'elec')&&(~isempty(MAL.elec)), set(handles.elec,'Value',MAL.elec); end
    
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
function varargout = log_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = cauto_Callback(h, eventdata, handles, varargin)
cauto=get(handles.cauto,'Value');
if cauto==1,
   set(handles.cmin,'Enable','inactive')
   set(handles.cmax,'Enable','inactive')
else
   set(handles.cmin,'Enable','on')
   set(handles.cmax,'Enable','on')
end

% --------------------------------------------------------------------
function varargout = cmin_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = cmax_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = cmap_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = ok_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.ok.
global MAL
MAL.cauto=get(handles.cauto,'Value');
MAL.clog=get(handles.log,'Value');
MAL.creverse=get(handles.creverse,'Value');
MAL.cmap=get(handles.cmap,'Value')-1;
MAL.cmax=str2num(get(handles.cmax,'String'));
MAL.cmin=str2num(get(handles.cmin,'String'));
MAL.style=get(handles.style,'Value')-1;
MAL.showgrid=get(handles.showgrid,'Value')-1;
MAL.high=get(handles.high,'Value')-1;
MAL.xdir=get(handles.xdir,'Value');
MAL.elec=get(handles.elec,'Value');
MAL.alpha=get(handles.alpha,'Value');
delete(gcbf);

% --------------------------------------------------------------------
function varargout = apply_Callback(h, eventdata, handles, varargin)
global MAL Mod
mal=MAL;
MAL.cauto=get(handles.cauto,'Value');
MAL.log=get(handles.log,'Value');
MAL.cmap=get(handles.cmap,'Value')-1;
MAL.cmax=str2num(get(handles.cmax,'String'));
MAL.cmin=str2num(get(handles.cmin,'String'));
MAL.style=get(handles.style,'Value')-1;
MAL.high=get(handles.high,'Value')-1;
MAL.xdir=get(handles.xdir,'Value');
MAL.elec=get(handles.elec,'Value');
MAL.alpha=get(handles.alpha,'Value');
%inv2d('showmod_Callback',gcbo,[],guidata(gcbo))
if isfield(MAL,'style')&&(MAL.style==3),
    set(figure(8),'MenuBar','none','NumberTitle','off','Name','LongProfileModel');
    plotlongmod(Mod,MAL);
else    
    draw2dmodel(Mod.x,Mod.z,Mod.M,MAL,Mod.Cov);
end

% --------------------------------------------------------------------
function varargout = cancel_Callback(h, eventdata, handles, varargin)
delete(gcbf);

% --------------------------------------------------------------------
function varargout = xdir_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = elec_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = cont_Callback(h, eventdata, handles, varargin)

% --- Executes during object creation, after setting all properties.
function style_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in style.
function style_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function high_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in high.
function high_Callback(hObject, eventdata, handles)


% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
ibrowse(['file://' pwd filesep 'doc' filesep 'options.html#graph'])


% --- Executes on button press in alpha.
function alpha_Callback(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of alpha


% --- Executes during object creation, after setting all properties.
function showgrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to showgrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in showgrid.
function showgrid_Callback(hObject, eventdata, handles)


% --- Executes on button press in creverse.
function creverse_Callback(hObject, eventdata, handles)
% hObject    handle to creverse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of creverse


