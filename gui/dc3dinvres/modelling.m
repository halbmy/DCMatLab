function varargout = modelling(varargin)
% MODELLING Application M-file for modelling.fig
%    FIG = MODELLING launch modelling GUI.
%    MODELLING('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 17-Nov-2006 11:55:12

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    iconify(fig);
    
    global Model
	set(handles.xmin,'String',num2str(min(Model.x)));
	set(handles.xmax,'String',num2str(max(Model.x)));
	set(handles.ymin,'String',num2str(min(Model.y)));
	set(handles.ymax,'String',num2str(max(Model.y)));
	set(handles.zmin,'String',num2str(min(Model.z)));
	set(handles.zmax,'String',num2str(max(Model.z)));
    if iscell(Model.M),
        return;
    else
        set(handles.rho,'String',num2str(round(Model.M(1))));
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
function varargout = xmin_Callback(h, eventdata, handles, varargin)
global Model
xx=str2num(get(h,'String'));
[xm,ind]=min(abs(Model.x-xx));
% ind=find(abs(Model.x-xx)==xm);
set(h,'String',num2str(Model.x(ind(1))));

% --------------------------------------------------------------------
function varargout = xmax_Callback(h, eventdata, handles, varargin)
global Model
xx=str2num(get(h,'String'));
[xm,ind]=min(abs(Model.x-xx));
% ind=find(abs(Model.x-xx)==xm);
set(h,'String',num2str(Model.x(ind(1))));

% --------------------------------------------------------------------
function varargout = ymin_Callback(h, eventdata, handles, varargin)
global Model
yy=str2num(get(h,'String'));
[ym,ind]=min(abs(Model.y-yy));
% ind=find(abs(Model.y-yy)==ym);
set(h,'String',num2str(Model.y(ind(1))));


% --------------------------------------------------------------------
function varargout = ymax_Callback(h, eventdata, handles, varargin)
global Model
yy=str2num(get(h,'String'));
[ym,ind]=min(abs(Model.y-yy));
% ind=find(abs(Model.y-yy)==ym);
set(h,'String',num2str(Model.y(ind(1))));

% --------------------------------------------------------------------
function varargout = zmin_Callback(h, eventdata, handles, varargin)
global Model
zz=str2num(get(h,'String'));
[zm,ind]=min(abs(Model.z-zz));
% ind=find(abs(Model.z-zz)==zm);
set(h,'String',num2str(Model.z(ind(1))));

% --------------------------------------------------------------------
function varargout = zmax_Callback(h, eventdata, handles, varargin)
global Model
zz=str2num(get(h,'String'));
[zm,ind]=min(abs(Model.z-zz));
% ind=find(abs(z-zz)==zm);
set(h,'String',num2str(Model.z(ind(1))));

% --------------------------------------------------------------------
function varargout = rho_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = all_Callback(h, eventdata, handles, varargin)
global Model
set(handles.xmin,'String',num2str(min(Model.x)));
set(handles.xmax,'String',num2str(max(Model.x)));
set(handles.ymin,'String',num2str(min(Model.y)));
set(handles.ymax,'String',num2str(max(Model.y)));
set(handles.zmin,'String',num2str(min(Model.z)));
set(handles.zmax,'String',num2str(max(Model.z)));


% --------------------------------------------------------------------
function varargout = set_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.all.
global Model
i1=find(Model.x==str2num(get(handles.xmin,'String')));
i2=find(Model.x==str2num(get(handles.xmax,'String')));
j1=find(Model.y==str2num(get(handles.ymin,'String')));
j2=find(Model.y==str2num(get(handles.ymax,'String')));
k1=find(Model.z==str2num(get(handles.zmin,'String')));
k2=find(Model.z==str2num(get(handles.zmax,'String')));
% if i1==1, i1=1; end
% if j1==1, j1=1; end
% if i2==length(x), i2=length(x); end
% if j2==length(y), j2=length(y); end
% if k2==length(z), k2=length(z); end

%fprintf('i=%d..%d j=%d..%d k=%d..%d\n',i1,i2,j1,j2,k1,k2);
rho=str2num(get(handles.rho,'String'));
ijk=(i2-i1)*(j2-j1)*(k2-k1);
if ijk>0,
    if rho(1)==0, set(handles.status,'String','Rho=0!');return; end
	Model.M(i1:i2-1,j1:j2-1,k1:k2-1)=rho(1);
	if (i1==0)&(j1==0),
        if k1==1, sq=1/rho(1); end
        Model.Bg(k1:k2)=rho;        
	end
	set(handles.status,'String',sprintf('Set %dx%dx%d=%d Cells to rho=%.1f',i2-i1,j2-j1,k2-k1,ijk,rho(1)));
else
    set(handles.status,'String',sprintf('No cells selected i=%d..%d j=%d..%d k=%d..%d',i1,i2,j1,j2,k1,k2));
end

% --------------------------------------------------------------------
function varargout = forward_Callback(h, eventdata, handles, varargin)
global N Model FOR
set(handles.status,'String','Solving forward problem with FD');
Model.R=mfdfwd3d(Model,N,FOR);
set(handles.status,'String',sprintf('Solved forward problem (RMS=%.2fperc.)',rms(N.r,Model.R)));

% --------------------------------------------------------------------
function varargout = show_Callback(h, eventdata, handles, varargin)
global Model MAL
set(figure(1),'Name','Model');
draw3dmodel(Model,MAL);

% --------------------------------------------------------------------
function varargout = close_Callback(h, eventdata, handles, varargin)
delete(gcbf);


% --- Executes on button press in fix.
function fix_Callback(hObject, eventdata, handles)
global Model FIX
i1=find(Model.x==str2num(get(handles.xmin,'String')));
i2=find(Model.x==str2num(get(handles.xmax,'String')));
j1=find(Model.y==str2num(get(handles.ymin,'String')));
j2=find(Model.y==str2num(get(handles.ymax,'String')));
k1=find(Model.z==str2num(get(handles.zmin,'String')));
k2=find(Model.z==str2num(get(handles.zmax,'String')));
ijk=(i2-i1)*(j2-j1)*(k2-k1);
if ijk>0,
    if ~isequal(size(Model.M),size(FIX)), FIX=zeros(size(Model.M)); end
	FIX(i1:i2-1,j1:j2-1,k1:k2-1)=-1;
	set(handles.status,'String',sprintf('Fixed %dx%dx%d=%d Cells',i2-i1,j2-j1,k2-k1,ijk));
else
    set(handles.status,'String',sprintf('No cells selected i=%d..%d j=%d..%d k=%d..%d',i1,i2,j1,j2,k1,k2));
end


% --- Executes on button press in compound.
function compound_Callback(hObject, eventdata, handles)
global Model FIX
i1=find(Model.x==str2num(get(handles.xmin,'String')));
i2=find(Model.x==str2num(get(handles.xmax,'String')));
j1=find(Model.y==str2num(get(handles.ymin,'String')));
j2=find(Model.y==str2num(get(handles.ymax,'String')));
k1=find(Model.z==str2num(get(handles.zmin,'String')));
k2=find(Model.z==str2num(get(handles.zmax,'String')));
ijk=(i2-i1)*(j2-j1)*(k2-k1);
if ijk>0,
    if ~isequal(size(Model.M),size(FIX)), FIX=zeros(size(Model.M)); end
	FIX(i1:i2-1,j1:j2-1,k1:k2-1)=max(FIX(:))+1;
	set(handles.status,'String',sprintf('Compounded %dx%dx%d=%d Cells',i2-i1,j2-j1,k2-k1,ijk));
else
    set(handles.status,'String',sprintf('No cells selected i=%d..%d j=%d..%d k=%d..%d',i1,i2,j1,j2,k1,k2));
end


% --- Executes on button press in decx.
function decx_Callback(hObject, eventdata, handles)
global Model XX
i1=find(Model.x==str2num(get(handles.xmin,'String')));
i2=find(Model.x==str2num(get(handles.xmax,'String')));
j1=find(Model.y==str2num(get(handles.ymin,'String')));
j2=find(Model.y==str2num(get(handles.ymax,'String')));
k1=find(Model.z==str2num(get(handles.zmin,'String')));
k2=find(Model.z==str2num(get(handles.zmax,'String')));
if i1==1, i1=2; end
if i2==length(Model.x), x2=length(Model.x)-1; end
ijk=(i2-i1+1)*(j2-j1)*(k2-k1);
if ijk>0,
    si=size(Model.M)-[1 0 0];
    if ~isequal(si,size(XX)), XX=zeros(si); end
	XX(i1-1:i2-1,j1:j2-1,k1:k2-1)=1;
	set(handles.status,'String',sprintf('Decoupled(x) %dx%dx%d=%d Faces',i2-i1+1,j2-j1,k2-k1,ijk));
else
    set(handles.status,'String',sprintf('No x-faces selected i=%d..%d j=%d..%d k=%d..%d',i1-1,i2-1,j1,j2,k1,k2));
end


% --- Executes on button press in decy.
function decy_Callback(hObject, eventdata, handles)
global Model YY
i1=find(Model.x==str2num(get(handles.xmin,'String')));
i2=find(Model.x==str2num(get(handles.xmax,'String')));
j1=find(Model.y==str2num(get(handles.ymin,'String')));
j2=find(Model.y==str2num(get(handles.ymax,'String')));
k1=find(Model.z==str2num(get(handles.zmin,'String')));
k2=find(Model.z==str2num(get(handles.zmax,'String')));
if j1==1, j1=2; end
if j2==length(Model.y), j2=length(Model.y)-1; end
ijk=(i2-i1)*(j2-j1+1)*(k2-k1);
if ijk>0,
    si=size(Model.M)-[0 1 0];
    if ~isequal(si,size(YY)), YY=zeros(si); end
	YY(i1:i2-1,j1-1:j2-1,k1:k2-1)=1;
	set(handles.status,'String',sprintf('Decoupled(y) %dx%dx%d=%d Faces',i2-i1,j2-j1+1,k2-k1,ijk));
else
    set(handles.status,'String',sprintf('No y-faces selected i=%d..%d j=%d..%d k=%d..%d',i1,i2,j1-1,j2-1,k1,k2));
end


% --- Executes on button press in decz.
function decz_Callback(hObject, eventdata, handles)
global Model ZZ
i1=find(Model.x==str2num(get(handles.xmin,'String')));
i2=find(Model.x==str2num(get(handles.xmax,'String')));
j1=find(Model.y==str2num(get(handles.ymin,'String')));
j2=find(Model.y==str2num(get(handles.ymax,'String')));
k1=find(Model.z==str2num(get(handles.zmin,'String')));
k2=find(Model.z==str2num(get(handles.zmax,'String')));
if k1<=1, k1=2; end
if k2>=length(Model.z), k2=length(Model.z)-1; end
ijk=(i2-i1)*(j2-j1)*(k2-k1+1);
if ijk>0,
    si=size(Model.M)-[0 0 1];
    if ~isequal(si,size(ZZ)), ZZ=zeros(si); end
	ZZ(i1:i2-1,j1:j2-1,k1-1:k2-1)=1;
	set(handles.status,'String',sprintf('Decoupled(z) %dx%dx%d=%d Faces',i2-i1,j2-j1,k2-k1+1,ijk));
else
    set(handles.status,'String',sprintf('No z-faces selected i=%d..%d j=%d..%d k=%d..%d',i1,i2,j1,j2,k1-1,k2-1));
end


% --- Executes on button press in removefix.
function removefix_Callback(hObject, eventdata, handles)
global FIX
FIX(:)=0;

% --- Executes on button press in removedec.
function removedec_Callback(hObject, eventdata, handles)
global XX YY ZZ
XX(:)=0;YY(:)=0;ZZ(:)=0;