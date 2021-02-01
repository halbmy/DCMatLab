function varargout = i_choose(varargin)
% I_CHOOSE Application M-file for i_choose.fig
%    FIG = I_SVD launch i_choose GUI.
%    I_CHOOSE('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 02-May-2004 11:24:09

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
    maxq=length(DM(1,:));
    load('rhoeta.mat');
    if(exist('ak')==1), qq=lkurv(rho,eta,ak); else qq=lkurv(rho,eta); end
    nl=0.1;
    if exist('ak')==1, nl=ak(qq); end
    set(handles.scut,'Max',maxq);
    set(handles.scut,'Min',minq);
    set(handles.scut,'Value',qq);
    set(handles.scut,'SliderStep',[1 10]/(maxq-minq));
    set(handles.scuttext,'String',num2str(qq));
    set(handles.regpar,'String',num2str(nl));
    
    i_choose('scut_Callback',0,[],handles);
%     malmodel(qq);

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
% if INV.redu==1
%     if INV.mitschicht==1
%         dM=[P zeros(size(P,1),size(PP,1));PP]'*dM;
%     else
%         dM=P'*dM;
%     end
% end
lbound=0;
if isfield(INV,'lbound')&&(INV.lbound>0), lbound=INV.lbound(1); end
if INV.lolo==1
    nM=(Mod.M(:)-lbound).*exp(dM(1:length(Mod.M(:))))+lbound;
else
    nM=Mod.M(:)+dM(1:length(Mod.M(:)));
end
draw2dmodel(Mod.x,Mod.z,nM,MAL);
    
% --------------------------------------------------------------------
function varargout = scut_Callback(h, eventdata, handles, varargin)
qq=round(get(handles.scut,'Value'));
set(handles.scut,'Value',qq);
set(handles.scuttext,'String',num2str(qq));
load('rhoeta.mat');
if exist('ak'), set(handles.regpar,'String',ak(qq)); end
malmodel(qq);
figure(handles.figure1);

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
if(exist('ak')), 
    warten=lkurv(rho,eta,ak);
else
    warten=lkurv(rho,eta);
end
delete(gcbf);


% --------------------------------------------------------------------
function varargout = regpar_Callback(h, eventdata, handles, varargin)
% global s U V ak
% nl=str2num(get(handles.regpar,'String'));
% load('rhoeta.mat');
% if exist('ak'),
%   qq=max(find(ak>nl));
%   if (~isempty(qq))&(qq>0),
%       set(handles.scut,'Value',qq);
%       set(handles.scuttext,'String',num2str(qq));
%       malmodel(qq);
%       figure(h);
%   end
% end
