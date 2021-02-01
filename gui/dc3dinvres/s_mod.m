function varargout = s_mod(varargin)
% S_MOD Application M-file for s_mod.fig
%    FIG = S_MOD launch s_mod GUI.
%    S_MOD('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 03-Nov-2005 12:25:49

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    iconify(fig);
    
    global Model N
    if isfield(Model,'x'), x=Model.x; else return; end
    if isfield(Model,'y'), y=Model.y; else return; end
    if isfield(Model,'z'), z=Model.z; else return; end
    if isfield(Model,'Bg'), Bg=Model.Bg; else return; end
    set(handles.editx,'String',sprintf('%g ',x));
    set(handles.edity,'String',sprintf('%g ',y));
    set(handles.editz,'String',sprintf('%g ',z));
    set(handles.layers,'String',num2str(length(z)-1));
    set(handles.zmax,'String',num2str(z(end)));
    set(handles.xmin,'String',num2str(min(x)));
    set(handles.xmax,'String',num2str(max(x)));
    set(handles.ymin,'String',num2str(min(y)));
    set(handles.ymax,'String',num2str(max(y)));

%     aa=sqrt(sum(diff(N.elec(:,1:2))'.^2));
%     D=min(aa(find(aa)));
%     if D<1, D=1; end
    dx=min(diff(x));if dx==0, dx=1; end
    dy=min(diff(y));if dy==0, dy=1; end
    set(handles.dx,'String',num2str(dx));
    set(handles.dy,'String',num2str(dy));
    set(handles.editbg,'String',sprintf('%d ',round(Bg)));
	% Wait for callbacks to run and window to be dismissed:
	uiwait(fig);

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
function varargout = editx_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit
X=str2num(get(handles.editx,'String'));
% set(handles.editx,'String',sprintf('%.1f ',X));


% --------------------------------------------------------------------
function varargout = edity_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol ha
Y=str2num(get(handles.edity,'String'));
set(handles.edity,'String',sprintf('%.1f ',Y));


% --------------------------------------------------------------------
function varargout = editz_Callback(h, eventdata, handles, varargin)
global Model
Z=str2num(get(handles.editz,'String'));
lbg=length(Z);
set(handles.editz,'String',sprintf('%g ',Z));
nBg=str2num(get(handles.editbg,'String'));
if length(nBg)==0, nBg=Model.Bg; end
if length(nBg)>lbg, nBg=nBg(1:lbg); end
while length(nBg)<lbg, nBg(end+1)=nBg(end); end
set(handles.editbg,'String',sprintf('%d ',round(nBg)));


% --------------------------------------------------------------------
function varargout = editbg_Callback(h, eventdata, handles, varargin)
global Model
Z=str2num(get(handles.editz,'String'));
lbg=length(Z);
nBg=str2num(get(handles.editbg,'String'));
if length(nBg)==0, nBg=Model.Bg; end
if length(nBg)>lbg, nBg=nBg(1:length(Model.Bg)); end
while length(nBg)<lbg, nBg(end+1)=nBg(end); end
set(handles.editbg,'String',sprintf('%d ',round(nBg)));


% --------------------------------------------------------------------
function varargout = dx_Callback(h, eventdata, handles, varargin)
global N
D=str2num(get(handles.dx,'String'));
minx=min(N.elec(:,1))-D;
maxx=max(N.elec(:,1))+D;
set(handles.xmin,'String',num2str(minx));
set(handles.xmax,'String',num2str(maxx));

% --------------------------------------------------------------------
function varargout = dy_Callback(h, eventdata, handles, varargin)
global N
D=str2num(get(handles.dy,'String'));
miny=min(N.elec(:,2))-D;
maxy=max(N.elec(:,2))+D;
set(handles.ymin,'String',num2str(miny));
set(handles.ymax,'String',num2str(maxy));

% --------------------------------------------------------------------
function varargout = ok_Callback(h, eventdata, handles, varargin)
global Model MAL
Model.x=str2num(get(handles.editx,'String'));
Model.y=str2num(get(handles.edity,'String'));
Model.z=str2num(get(handles.editz,'String'));
Model.Bg=str2num(get(handles.editbg,'String'));
I=length(Model.x)-1;J=length(Model.y)-1;K=length(Model.z)-1;
if ~isequal(size(Model.M),[I J K]),
    Model.M=Model.Bg(1)*ones(I,J,K);
end
if ~isequal(Model.Bg,0),
    for l = 1:length(Model.Bg)-1
        Model.M(:,:,l)=Model.Bg(l);
    end
end
delete(gcbf);

% --------------------------------------------------------------------
function varargout = cancel_Callback(h, eventdata, handles, varargin)
delete(gcbf);

% --------------------------------------------------------------------
function varargout = layers_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = zpar_Callback(h, eventdata, handles, varargin)
global N
lay=round(str2num(get(handles.layers,'String')));
if lay>1,
    z=zparam(N,lay);
    Bg=str2num(get(handles.editbg,'String'));
    if length(Bg)>=length(z)-1, Bg=Bg(1:length(z)-1); end
    if length(Bg)<length(z)-1, Bg(end:length(z)-1)=Bg(end); end
    set(handles.editz,'String',sprintf('%.1f ',z));
    set(handles.editbg,'String',sprintf('%d ',round(Bg)));
end

% --------------------------------------------------------------------
function varargout = linear_Callback(h, eventdata, handles, varargin)


% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
ibrowse(['file://' pwd filesep 'doc' filesep 'model.html#modpar']) 

% --- Executes during object creation, after setting all properties.
function zmax_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function zmax_Callback(hObject, eventdata, handles)

% --- Executes on button press in log.
function log_Callback(hObject, eventdata, handles)
global N
lay=round(str2num(get(handles.layers,'String')));
zmax=str2num(get(handles.zmax,'String'));
Bg=str2num(get(handles.editbg,'String'));
if lay>1,
    z=zparam(N,lay,zmax);
    if length(Bg)>length(z), Bg=Bg(1:length(z)); end
    if length(Bg)<length(z), Bg(end:length(z))=Bg(end); end
    set(handles.editz,'String',sprintf('%g ',z));
    set(handles.editbg,'String',sprintf('%d ',round(Bg)));
end

% --- Executes on button press in xmake.
function xmake_Callback(hObject, eventdata, handles)
xmin=str2num(get(handles.xmin,'String'));
xmax=str2num(get(handles.xmax,'String'));
dx=str2num(get(handles.dx,'String'));
x=xmin:dx:xmax;
set(handles.editx,'String',sprintf('%g ',x));

% --- Executes on button press in ymake.
function ymake_Callback(hObject, eventdata, handles)
ymin=str2num(get(handles.ymin,'String'));
ymax=str2num(get(handles.ymax,'String'));
dy=str2num(get(handles.dy,'String'));
y=ymin:dy:ymax;
set(handles.edity,'String',sprintf('%g ',y));

% --- Executes during object creation, after setting all properties.
function xmin_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function xmin_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function xmax_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function xmax_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ymin_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function ymin_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ymax_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function ymax_Callback(hObject, eventdata, handles)


% --- Executes on button press in lin.
function lin_Callback(hObject, eventdata, handles)
lay=round(str2num(get(handles.layers,'String')));
zmax=str2num(get(handles.zmax,'String'));
Bg=str2num(get(handles.editbg,'String'));
if lay>1,
    z=round(linspace(0,zmax,lay+1)*1000)/1000;
    if length(Bg)>length(z), Bg=Bg(1:length(z)); end
    if length(Bg)<length(z), Bg(end:length(z))=Bg(end); end
    set(handles.editz,'String',sprintf('%g ',z));
    set(handles.editbg,'String',sprintf('%d ',round(Bg)));
end
