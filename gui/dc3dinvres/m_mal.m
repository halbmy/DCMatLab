function varargout = s_mal(varargin)
% M_MAL Application M-file for s_mal.fig
%    FIG = S_MAL launch s_mal GUI.
%    S_MAL('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 16-Feb-2004 13:40:14

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    iconify(fig);
   
    global MAL
    save('s_mal.mat','MAL');
    
    if isfield(MAL,'vie'), set(handles.vie,'Value',MAL.vie+1); end
    if isfield(MAL,'xy'), set(handles.xy,'Value',MAL.xy); end
    if isfield(MAL,'xdir'), set(handles.xdir,'Value',MAL.xdir); end
    if isfield(MAL,'ydir'), set(handles.ydir,'Value',MAL.ydir); end
    if isfield(MAL,'log'), set(handles.log,'Value',MAL.log); end
    if isfield(MAL,'clog'), set(handles.log,'Value',MAL.clog); end
    if isfield(MAL,'cauto'), set(handles.cauto,'Value',MAL.cauto); end
    if isfield(MAL,'cmap'), set(handles.cmap,'Value',MAL.cmap); end
    if isfield(MAL,'cmin'), set(handles.cmin,'String',num2str(MAL.cmin)); end
    if isfield(MAL,'cmax'), set(handles.cmax,'String',num2str(MAL.cmax)); end
    if isfield(MAL,'cauto')&&(MAL.cauto==1),
        set(handles.cmin,'Enable','inactive')
        set(handles.cmax,'Enable','inactive')
    end
    if isfield(MAL,'nu'), set(handles.nu,'String',num2str(MAL.nu)); end
    if isfield(MAL,'nv'), set(handles.nv,'String',num2str(MAL.nv)); end
    if isfield(MAL,'startwith'), set(handles.startwith,'String',num2str(MAL.startwith)); end
    if isfield(MAL,'elec'), set(handles.elec,'String',num2str(MAL.elec)); end
    if isfield(MAL,'cont'), set(handles.cont,'String',strrep(num2str(MAL.cont),'  ',' ')); end
    if isfield(MAL,'alpha'), set(handles.alpha,'Value',MAL.alpha>0); end
    
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
function varargout = xy_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.xy.
global MAL
MAL.xy=get(handles.xy,'Value');

% --------------------------------------------------------------------
function varargout = log_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.log.
global MAL
MAL.log=get(handles.log,'Value');
MAL.clog=get(handles.log,'Value');

% --------------------------------------------------------------------
function varargout = cauto_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.log.
global MAL
MAL.cauto=get(handles.cauto,'Value');
if MAL.cauto==1,
   set(handles.cmin,'Enable','inactive')
   set(handles.cmax,'Enable','inactive')
else
   set(handles.cmin,'Enable','on')
   set(handles.cmax,'Enable','on')
end

% --------------------------------------------------------------------
function varargout = cmin_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.cmin.
global MAL
MAL.cmin=str2num(get(handles.cmin,'String'));

% --------------------------------------------------------------------
function varargout = cmax_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.cmax.
global MAL
MAL.cmax=str2num(get(handles.cmax,'String'));

% --------------------------------------------------------------------
function varargout = cmap_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.popupmenu1.
global MAL
MAL.cmap=get(handles.cmap,'Value');

% --------------------------------------------------------------------
function varargout = nu_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.nu.
global MAL
MAL.nu=str2num(get(handles.nu,'String'));

% --------------------------------------------------------------------
function varargout = nv_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.nv.
global MAL
MAL.nv=str2num(get(handles.nv,'String'));

% --------------------------------------------------------------------
function varargout = startwith_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit5.
global MAL
MAL.startwith=str2num(get(handles.startwith,'String'));

% --------------------------------------------------------------------
function varargout = ok_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.ok.
global MAL
delete('s_mal.mat');
delete(gcbf);

% --------------------------------------------------------------------
function varargout = apply_Callback(h, eventdata, handles, varargin)
global Model MAL Cov
dc3dinvres('showmod_Callback',gcbo,[],guidata(gcbo))
% set(figure(1),'NumberTitle','off','Name','Model');
% draw3dmodel(Model,MAL,[],Cov);

% --------------------------------------------------------------------
function varargout = cancel_Callback(h, eventdata, handles, varargin)
global MAL
load('s_mal.mat');
delete('s_mal.mat');
delete(gcbf);

% --------------------------------------------------------------------
function varargout = xdir_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.xdir.
global MAL
MAL.xdir=get(handles.xdir,'Value');

% --------------------------------------------------------------------
function varargout = ydir_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.ydir.
global MAL
MAL.ydir=get(handles.ydir,'Value');

% --------------------------------------------------------------------
function varargout = elec_Callback(h, eventdata, handles, varargin)
global MAL
MAL.elec=str2num(get(handles.elec,'String'));

% --------------------------------------------------------------------
function varargout = cont_Callback(h, eventdata, handles, varargin)
global MAL
MAL.cont=str2num(get(handles.cont,'String'));

% --------------------------------------------------------------------
function varargout = stdcont_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.stdcont.
global MAL
MAL.cont=[10 20 50 100 200 500 1000];
set(handles.cont,'String',strrep(num2str(MAL.cont),'  ',' '));


% --------------------------------------------------------------------
function varargout = vie_Callback(h, eventdata, handles, varargin)
global MAL
MAL.vie=get(handles.vie,'Value')-1;


% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
ibrowse(['file://' pwd filesep '3ddoc' filesep 'options.html#graph']) 


% --- Executes on button press in alpha.
function alpha_Callback(hObject, eventdata, handles)
global MAL
MAL.alpha=get(handles.alpha,'Value');

