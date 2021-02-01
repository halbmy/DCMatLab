function varargout = s_inv(varargin)
% S_INV Application M-file for s_inv.fig
%    FIG = S_INV launch s_inv GUI.
%    S_INV('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 10-Sep-2003 18:02:06

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    iconify(fig);
    ss=get(handles.method,'String');
    ss{3}='Reverse';
    set(handles.method,'String',ss);
    global FOR
    if isfield(FOR,'method'), set(handles.method,'Value',FOR.method+1); end
    if isfield(FOR,'refine'), set(handles.refining,'Value',min(FOR.refine+1,5)); end
    if isfield(FOR,'zref'), set(handles.zref,'Value',FOR.zref); end
    if isfield(FOR,'rand'), set(handles.rand,'String',num2str(FOR.rand)); end
    if isfield(FOR,'zusatz'), set(handles.zusatz,'String',num2str(FOR.zusatz)); end
    if isfield(FOR,'prolong'), set(handles.prolong,'String',num2str(FOR.prolong)); end
    if isfield(FOR,'fillup'), set(handles.fillup,'Value',FOR.fillup+1); end
    %if isfield(FOR,'acc'), set(handles.acc,'String',num2str(FOR.acc)); end
    %if isfield(FOR,'tol'), set(handles.tol,'String',num2str(FOR.tol)); end
    %if isfield(FOR,'maxit'), set(handles.maxit,'String',num2str(FOR.maxit)); end
	if nargout > 0, varargout{1} = fig; end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end

end

% --------------------------------------------------------------------
function varargout = inv_ok_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.inv_ok.
global FOR
% mmm=str2num(get(handles.acc,'String'));
% if mmm>0, FOR.acc=mmm; end
% mmm=str2num(get(handles.tol,'String'));
% if mmm>0, FOR.tol=mmm; end
% mmm=fix(str2num(get(handles.maxit,'String')));
% if mmm>0, FOR.maxit=mmm; end
 mmm=fix(str2num(get(handles.rand,'String')));
if mmm>0, FOR.rand=mmm; end
mmm=fix(str2num(get(handles.prolong,'String')));
if mmm>0, FOR.prolong=mmm; end
mmm=fix(str2num(get(handles.zusatz,'String')));
if mmm>0, FOR.zusatz=mmm; end
FOR.method=get(handles.method,'Value')-1;
FOR.refine=get(handles.refining,'Value')-1;
if FOR.refine==4, FOR.refine=5; end
FOR.fillup=get(handles.fillup,'Value')-1;
FOR.zref=get(handles.zref,'Value');
delete(gcbf);

% --------------------------------------------------------------------
function varargout = inv_cancel_Callback(h, eventdata, handles, varargin)
delete(gcbf);

% --------------------------------------------------------------------
function varargout = refining_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = calc_Callback(h, eventdata, handles, varargin)
global Mod N FOR INV RMS CHIQ
newfor=FOR;
% mmm=str2num(get(handles.acc,'String'));
% if mmm>0, newfor.acc=mmm; end
% mmm=str2num(get(handles.tol,'String'));
% if mmm>0, newfor.tol=mmm; end
% mmm=fix(str2num(get(handles.maxit,'String')));
% if mmm>0, newfor.maxit=mmm; end
mmm=fix(str2num(get(handles.rand,'String')));
if mmm>0, newfor.rand=mmm; end
mmm=fix(str2num(get(handles.prolong,'String')));
if mmm>0, new.prolong=mmm; end
mmm=fix(str2num(get(handles.zusatz,'String')));
if mmm>0, newfor.zusatz=mmm; end
newfor.method=get(handles.method,'Value')-1;
newfor.refine=get(handles.refining,'Value')-1;
newfor.zref=get(handles.zref,'Value');
newfor.fillup=get(handles.fillup,'Value')-1;

Mod.R=abs(dcfwd2d(Mod.x,Mod.z,Mod.M,Mod.Lay,N,newfor));
RMS=[RMS rms(N.r,Mod.R,INV.lolo)];
CHIQ=[CHIQ chi2(N.r,Mod.R,N.err)];
message(['CHI^2 = ',sprintf('%.2f ',CHIQ) '(RMS=' sprintf('%.2f',RMS(end)) ')']);
% message(strcat('RMS = ',sprintf('%.2f ',RMS)));


% --- Executes during object creation, after setting all properties.
function fillup_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in fillup.
function fillup_Callback(hObject, eventdata, handles)


% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
ibrowse(['file://' pwd filesep 'doc' filesep 'options.html#for'])


