function varargout = choose2dmodel(varargin)
% CHOOSE2DMODEL Application Mod.M-file for choose2dmodel.fig
%    FIG = I_SVD launch choose2dmodel GUI.
%    CHOOSE2DMODEL('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 06-Jan-2004 14:59:23

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    iconify(fig);

	% Wait for callbacks to run and window to be dismissed:

	if nargout > 0
		varargout{1} = fig;
    end
    
    global Mod P dM DM
    minq=1;
    maxq=size(DM,2);
    load('rhoeta.mat');
    if(exist('ak','var')==1), 
        qq=lkurv(rho,eta,ak); 
    else 
        fprintf('no lambdas present\n');
        qq=lkurv(rho,eta); 
    end
    nl=1;
    if exist('ak','var')==1, nl=ak(qq); end
    set(handles.scut,'Max',maxq);
    set(handles.scut,'Min',minq);
    set(handles.scut,'Value',qq);
    set(handles.scut,'SliderStep',[1 10]/(maxq-minq));
    set(handles.scuttext,'String',num2str(qq));
    set(handles.lambda,'String',num2str(nl));
    malmodel(qq);

    uiwait(fig);
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		 [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		%tom = feval(varargin{:}); % FEVAL switchyard
		%if tom==1, [varargout{1:nargout}] = tom; end
	catch
		disp(lasterr);
	end
end

function malmodel(qq)
global DM Mod INV MAL P PP
dM=DM(:,qq);
if INV.redu==1
    if INV.mitschicht==1
        dM=[P zeros(size(P,1),size(PP,1));PP]'*dM;
    else
        dM=P'*dM;
    end
end
%if INV.lolo==1
nM=Mod.M(:).*exp(dM(1:length(Mod.M(:))));
%else
%    nM=Mod.M(:)+dM(1:length(Mod.M(:)));
%end
draw2dmodel(Mod.x,Mod.z,nM,MAL);
    
% --------------------------------------------------------------------
function varargout = scut_Callback(h, eventdata, handles, varargin)
qq=round(get(handles.scut,'Value'));
set(handles.scut,'Value',qq);
set(handles.scuttext,'String',num2str(qq));
load('rhoeta.mat');
if exist('ak')==1, 
    set(handles.lambda,'String',num2str(ak(qq))); 
else
    set(handles.lambda,'String',num2str(qq)); 
end
malmodel(qq);

% --------------------------------------------------------------------
function varargout = ok_Callback(h, eventdata, handles, varargin)
global warten Mod
warten=round(get(handles.scut,'Value'));
delete(gcbf);

% --------------------------------------------------------------------
function varargout = video_Callback(h, eventdata, handles, varargin)
maxq=get(handles.scut,'Max');
figure;
Movie=moviein(maxq);
for l = 1:maxq,
    malmodel(l);
    drawnow;
    Movie(:,l)=getframe(gcf);
end
%movie2avi(Movie,'movie.avi','fps',1,'keyframe',1,'quality',100);

% --------------------------------------------------------------------
function varargout = cancel_Callback(h, eventdata, handles, varargin)
global Mod warten
load('rhoeta.mat');
if exist('ak')==1, 
    warten=lkurv(rho,eta,ak);
else
    warten=lkurv(rho,eta);
end
delete(gcbf);

% --------------------------------------------------------------------
function varargout = noiselevel_Callback(h, eventdata, handles, varargin)
nl=str2double(get(handles.noiselevel,'String'));
load('rhoeta.mat');
if exist('ak','var')==1,
  qq=max(find(ak>nl));
  if (~isempty(qq))&(qq>0),
      set(handles.scut,'Value',qq);
      set(handles.scuttext,'String',num2str(qq));
      malmodel(qq);
      figure(h);
  end
end

% --- Executes during object creation, after setting all properties.
function lambda_CreateFcn(hObject, eventdata, handles)

function lambda_Callback(hObject, eventdata, handles)
nl=str2double(get(handles.lambda,'String'));
load('rhoeta.mat');
if exist('ak','var')==1,
  qq=max(find(ak>nl));
  if (~isempty(qq))&&(qq>0),
      set(handles.scut,'Value',qq);
      set(handles.scuttext,'String',num2str(qq));
      malmodel(qq);
  end
end
