function varargout = s_mod(varargin)
% S_MOD Application Mod.M-file for s_mod.fig
%    FIG = S_MOD launch s_mod GUI.
%    S_MOD('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 28-Nov-2005 17:53:15

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    iconify(fig);
    
    global Mod N
    set(handles.editx,'String',sprintf('%g ',Mod.x));
    set(handles.editz,'String',sprintf('%g ',Mod.z));
    set(handles.layers,'String',num2str(length(Mod.z)-1));
    set(handles.zmax,'String',num2str(max(Mod.z)));
    
    aa=sqrt(sum(diff(N.elec(:,1:2))'.^2));
    D=round(min(aa(find(aa)))*20)/20;%5cm
    if D<0.1, D=1; end

    set(handles.dx,'String',num2str(D));
    set(handles.minx,'String',num2str(min(Mod.x)));
    set(handles.maxx,'String',num2str(max(Mod.x)));
    set(handles.editbg,'String',sprintf('%g ',rndig(Mod.Lay)));
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
set(handles.editx,'String',sprintf('%g ',X));



% --------------------------------------------------------------------
function varargout = editz_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.editz.
Z=str2num(get(handles.editz,'String'));
set(handles.editz,'String',sprintf('%g ',Z));

% --------------------------------------------------------------------
function varargout = editbg_Callback(h, eventdata, handles, varargin)
global Mod
Z=str2num(get(handles.editz,'String'));
lbg=length(Z);
nBg=str2num(get(handles.editbg,'String'));
if isempty(nBg), nBg=Mod.Lay; end
if length(nBg)>lbg, nBg=nBg(1:length(Mod.Lay)); end
while length(nBg)<lbg, nBg(end+1)=nBg(end); end
set(handles.editbg,'String',sprintf('%g ',rndig(nBg)));

% --------------------------------------------------------------------
function varargout = dx_Callback(h, eventdata, handles, varargin)
global N


% --------------------------------------------------------------------
function varargout = ok_Callback(h, eventdata, handles, varargin)

global Mod MAL S
Mod.x=str2num(get(handles.editx,'String'));
Mod.z=str2num(get(handles.editz,'String'));
Mod.Lay=str2num(get(handles.editbg,'String'));
I=length(Mod.x)-1;K=length(Mod.z)-1;
if Mod.Lay(1)==0, 
    Mod.M=Mod.M(1,1)*ones(I,K);
else
    Mod.M=Mod.Lay(1)*ones(I,K);
end
for l = 2:length(Mod.Lay)-1,
    if Mod.Lay(l)>0, Mod.M(:,l)=Mod.Lay(l); end
end
delete(gcbf);

% --------------------------------------------------------------------
function varargout = cancel_Callback(h, eventdata, handles, varargin)
delete(gcbf);

% --------------------------------------------------------------------
function varargout = zpar_Callback(h, eventdata, handles, varargin)
global N Mod
lay=round(str2num(get(handles.layers,'String')));
if lay>1,
    Mod.z=zparam(N,lay,0);
    if length(Mod.Lay)>=length(Mod.z)-1, Mod.Lay=Mod.Lay(1:length(Mod.z)-1); end
    if length(Mod.Lay)<length(Mod.z)-1, Mod.Lay(end:length(Mod.z)-1)=Mod.Lay(end); end
    set(handles.editz,'String',sprintf('%g ',Mod.z));
    set(handles.editbg,'String',sprintf('%g ',rndig(Mod.Lay)));
end

% --------------------------------------------------------------------
function varargout = logsp_Callback(h, eventdata, handles, varargin)
global N Mod
lay=round(str2num(get(handles.layers,'String')));
zmax=str2num(get(handles.zmax,'String'));
if lay>1,
    Mod.z=zparam(N,lay,zmax);
    if length(Mod.Lay)>=length(Mod.z)-1, Mod.Lay=Mod.Lay(1:length(Mod.z)-1); end
    if length(Mod.Lay)<length(Mod.z)-1, Mod.Lay(end:length(Mod.z)-1)=Mod.Lay(end); end
    set(handles.editz,'String',sprintf('%g ',Mod.z));
    set(handles.editbg,'String',sprintf('%g ',rndig(Mod.Lay)));
end

% --------------------------------------------------------------------
function varargout = linear_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = layers_Callback(h, eventdata, handles, varargin)


% --- Executes during object creation, after setting all properties.
function zmax_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function zmax_Callback(hObject, eventdata, handles)


% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
ibrowse(['file://' pwd filesep 'doc' filesep 'model.html#modpar'])


% --- Executes during object creation, after setting all properties.
function minx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function minx_Callback(hObject, eventdata, handles)
% hObject    handle to minx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minx as text
%        str2double(get(hObject,'String')) returns contents of minx as a double


% --- Executes during object creation, after setting all properties.
function maxx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function maxx_Callback(hObject, eventdata, handles)
% hObject    handle to maxx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxx as text
%        str2double(get(hObject,'String')) returns contents of maxx as a double


% --- Executes on button press in makex.
function makex_Callback(hObject, eventdata, handles)
D=str2num(get(handles.dx,'String'));
minx=str2num(get(handles.minx,'String'));
maxx=str2num(get(handles.maxx,'String'));
zusatz=0;
X=minx-zusatz*D:D:maxx+zusatz*D;
set(handles.editx,'String',sprintf('%g ',X));


% --- Executes on button press in linsp.
function linsp_Callback(hObject, eventdata, handles)
lay=round(str2num(get(handles.layers,'String')));
zmax=str2num(get(handles.zmax,'String'));
Mod.Lay=str2num(get(handles.editbg,'String'));
if lay>1,
    Mod.z=rndig(linspace(0,zmax,lay+1),3);
    if length(Mod.Lay)>length(Mod.z), Mod.Lay=Mod.Lay(1:length(Mod.z)); end
    if length(Mod.Lay)<length(Mod.z), Mod.Lay(end:length(Mod.z))=Mod.Lay(end); end
    set(handles.editz,'String',sprintf('%g ',Mod.z));
    set(handles.editbg,'String',sprintf('%g ',rndig(Mod.Lay)));
end 
