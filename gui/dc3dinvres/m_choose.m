function varargout = i_choose(varargin)
% I_CHOOSE Application M-file for i_choose.fig
%    FIG = I_SVD launch i_choose GUI.
%    I_CHOOSE('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 03-Sep-2001 12:59:46

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
    
    global Model P DM
    minq=1;
    maxq=length(DM(1,:));
    load rhoeta
    if exist('ak')==1, qq=lkurv(rho,eta,ak); else qq=lkurv(rho,eta); end
    nl=0.1;
    if exist('ak')==1, nl=ak(qq); end
    set(handles.scut,'Max',maxq);
    set(handles.scut,'Min',minq);
    set(handles.scut,'Value',qq);
    set(handles.scut,'SliderStep',[1 10]/(maxq-minq));
    set(handles.scuttext,'String',num2str(qq));
    set(handles.noiselevel,'String',num2str(nl));
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

global P Model DM MAL PP INV
dM=DM(:,qq);
if INV.redu==1
    if INV.mitschicht==1
        dM=[P zeros(length(P(:,1)),length(PP(:,1)));PP]'*dM;
    else
        dM=P'*dM;
    end
end
nM=modelupdate(Model,dM,2-(INV.lolo>0));;
figure(1);
draw3dmodel(nM,MAL);

% --------------------------------------------------------------------
function varargout = scut_Callback(h, eventdata, handles, varargin)

qq=round(get(handles.scut,'Value'));
set(handles.scut,'Value',qq);
set(handles.scuttext,'String',num2str(qq));
load rhoeta
if exist('ak')==1, set(handles.noiselevel,'String',ak(qq)); end
malmodel(qq);

% --------------------------------------------------------------------
function varargout = ok_Callback(h, eventdata, handles, varargin)

global warten
warten=round(get(handles.scut,'Value'));
delete(gcbf);


% --------------------------------------------------------------------
function varargout = video_Callback(h, eventdata, handles, varargin)
maxq=get(handles.scut,'Max');
% Movie=moviein(maxq);
% for l = 1:maxq,
%     malmodel(l);
%     Movie(:,l)=getframe;
%     makeframe(l);
% end
% makemovie([1:maxq],'choose',1,0,'gif');
% aviobj=avifile('choose1.avi','Quality',100,'Fps',1)
% for l = 1:maxq,
%     malmodel(l);
%     pause(1.0);
%     frame=getframe(gcf);
%     aviobj=addframe(aviobj,frame);
% end
% aviobj=close(aviobj);
%figure(1);
Movie=moviein(maxq);
for l = 1:maxq,
    malmodel(l);
    pause(0.5);
    Movie(:,l)=getframe(gcf);
end
%movie2avi(Movie,'movie.avi','fps',1,'keyframe',1,'quality',100);

% --------------------------------------------------------------------
function varargout = cancel_Callback(h, eventdata, handles, varargin)
global Model
load rhoeta
if exist('ak')==1, warten=lkurv(rho,eta,ak); else warten=lkurv(rho,eta); end
delete(gcbf);


% --------------------------------------------------------------------
function varargout = noiselevel_Callback(h, eventdata, handles, varargin)

global s U V

nl=str2num(get(handles.noiselevel,'String'));
load rhoeta
if exist('ak')==1,
    qq=max(find(ak>nl));
    if (~isempty(qq))&(qq>0),
        set(handles.scut,'Value',qq);
        set(handles.scuttext,'String',num2str(qq));
        malmodel(qq);
    end
end
