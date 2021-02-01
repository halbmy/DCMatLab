function varargout = s_for(varargin)
% S_FOR Application M-file for s_inv.fig
%    FIG = S_FOR launch s_for GUI.
%    S_FOR('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 25-May-2004 15:20:39

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    iconify(fig);
    global FOR
    %save('m_for.mat','FOR');
    if isfield(FOR,'refine'), set(handles.refining,'Value',FOR.refine+1); end
    if isfield(FOR,'zref'), set(handles.zref,'Value',FOR.zref); end
    if isfield(FOR,'rand'), set(handles.rand,'String',num2str(FOR.rand)); end
    if isfield(FOR,'zusatz'), set(handles.zusatz,'String',num2str(FOR.zusatz)); end
    if isfield(FOR,'prolong'), set(handles.prolong,'String',num2str(FOR.prolong)); end
    if isfield(FOR,'acc'), set(handles.acc,'String',num2str(FOR.acc)); end
    if isfield(FOR,'tol'), set(handles.tol,'String',num2str(FOR.tol)); end
    if isfield(FOR,'maxit'), set(handles.maxit,'String',num2str(FOR.maxit)); end
    if isfield(FOR,'method'), set(handles.method,'Value',FOR.method+1); end
    if isfield(FOR,'direct'), 
        set(handles.direct,'Value',FOR.direct+2); 
        if FOR.direct==1, 
            set(handles.acc,'Enable','off');
            set(handles.tol,'Enable','off');
            set(handles.maxit,'Enable','off');
        end
    end
    if isfield(FOR,'fillup'), set(handles.fillup,'Value',(FOR.fillup>0)+1); end
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
function varargout = ok_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.ok.
global FOR
mmm=str2num(get(handles.acc,'String'));
if mmm>0
    FOR.acc=mmm;
end
mmm=str2num(get(handles.tol,'String'));
if mmm>0
    FOR.tol=mmm;
end
mmm=fix(str2num(get(handles.maxit,'String')));
if mmm>0
    FOR.maxit=mmm;
end
mmm=fix(str2num(get(handles.rand,'String')));
if mmm>0
    FOR.rand=mmm;
end
mmm=fix(str2num(get(handles.prolong,'String')));
if mmm>0
    FOR.prolong=mmm;
end
mmm=fix(str2num(get(handles.zusatz,'String')));
if mmm>0
    FOR.zusatz=mmm;
end
FOR.method=get(handles.method,'Value')-1;
FOR.refine=get(handles.refining,'Value')-1;
FOR.zref=get(handles.zref,'Value');
FOR.direct=get(handles.direct,'Value')-2;
FOR.fillup=get(handles.fillup,'Value')-1;
delete(gcbf);

% --------------------------------------------------------------------
function varargout = cancel_Callback(h, eventdata, handles, varargin)
delete(gcbf);

% --------------------------------------------------------------------
function varargout = listbox1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.listbox1.
disp('listbox1 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = refining_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = direct_Callback(h, eventdata, handles, varargin)
onoff='on';
if get(handles.direct,'Value')>2, onoff='off'; end
set(handles.acc,'Enable',onoff);
set(handles.tol,'Enable',onoff);
set(handles.maxit,'Enable',onoff);

% --------------------------------------------------------------------
function varargout = calc_Callback(h, eventdata, handles, varargin)
global FOR
optfor=FOR;
mmm=str2num(get(handles.acc,'String'));
if mmm>0
    optfor.acc=mmm;
end
mmm=str2num(get(handles.tol,'String'));
if mmm>0
    optfor.tol=mmm;
end
mmm=fix(str2num(get(handles.maxit,'String')));
if mmm>0
    optfor.maxit=mmm;
end
mmm=fix(str2num(get(handles.rand,'String')));
if mmm>0
    optfor.rand=mmm;
end
mmm=fix(str2num(get(handles.prolong,'String')));
if mmm>0
    optfor.prolong=mmm;
end
mmm=fix(str2num(get(handles.zusatz,'String')));
if mmm>0
    optfor.zusatz=mmm;
end
optfor.method=get(handles.method,'Value')-1;
optfor.refine=get(handles.refining,'Value')-1;
optfor.zref=get(handles.zref,'Value');
optfor.direct=get(handles.direct,'Value');
optfor.fillup=get(handles.fillup,'Value')-1;

global Model N R RMS CHIQ INV
R=mfdfwd3d(Model,N,optfor);
RMS=[RMS rms(N.r,R,INV.lolo)];
CHIQ=[CHIQ chi2(N.r,R,N.err,INV.lolo)];
message(['CHI^2 = ',sprintf('%.2f ',CHIQ) '(RMS=' sprintf('%.2f',RMS(end)) ')']); 

% --- Executes during object creation, after setting all properties.
function zref_CreateFcn(hObject, eventdata, handles)
% if ispc
%     set(hObject,'BackgroundColor','white');
% else
%     set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
% end


% --- Executes on selection change in zref.
function zref_Callback(hObject, eventdata, handles)


% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
ibrowse(['file://' pwd filesep 'doc' filesep 'options.html#for']) 


% --- Executes during object creation, after setting all properties.
function fillup_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in fillup.
function fillup_Callback(hObject, eventdata, handles)

