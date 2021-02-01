function erg = m_svd(varargin)
% M_SVD Application M-file for m_svd.fig
%    FIG = M_SVD launch m_svd GUI.
%    M_SVD('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 13-Mar-2001 13:45:17

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    iconify(fig);

	%if nargout > 0
	%	varargout{1} = fig;
	%end
    
    global s VD VM Model P dM
    minq=1;
    nl=1*0.01;
    maxq=max(find(s>max(s)*nl));
    nl=10*0.01;
    qq=max(find(s>max(s)*nl));
    set(handles.scut,'Min',minq);
    set(handles.scut,'Max',maxq);
    set(handles.scut,'Value',qq);
    set(handles.scut,'SliderStep',[1 10]/(maxq-minq));
    set(handles.scuttext,'String',num2str(qq));
    set(handles.noiselevel,'String',num2str(nl));
    malmodel(qq);
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		 [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		%tom = feval(varargin{:}); % FEVAL switchyard
		%if tom==1, [varargout{1:nargout}] = tom; end
	catch
		disp(lasterr);
	end
end
    
% --------------------------------------------------------------------
function varargout = scut_Callback(h, eventdata, handles, varargin)

global s

qq=round(get(handles.scut,'Value'));
set(handles.scut,'Value',qq);
set(handles.scuttext,'String',num2str(qq));
set(handles.noiselevel,'String',fix(100*s(qq)/s(1))/100);
malmodel(qq);

% --------------------------------------------------------------------
function varargout = ok_Callback(h, eventdata, handles, varargin)

global warten
warten=round(get(handles.scut,'Value'));
delete(gcbf);

% --------------------------------------------------------------------
function varargout = video_Callback(h, eventdata, handles, varargin)
maxq=get(handles.scut,'Max');
Movie=moviein(maxq);
for l = 1:maxq,
    malmodel(l);
    Movie(:,l)=getframe(gcf);
end
%movie2avi(Movie,'movie.avi','fps',1,'keyframe',1,'quality',100);

% --------------------------------------------------------------------
function varargout = cancel_Callback(h, eventdata, handles, varargin)
global M warten s
nl=10*0.01;
warten=max(find(s>max(s)*nl));
delete(gcbf);


% --------------------------------------------------------------------
function varargout = noiselevel_Callback(h, eventdata, handles, varargin)
global s

nl=str2num(get(handles.noiselevel,'String'));
qq=max(find(s>(s(1)*nl(1))));
if ~isempty(qq)
    set(handles.scut,'Value',qq);
    set(handles.scuttext,'String',num2str(qq));
    malmodel(qq);
end

function malmodel(qq)

global Model VD VM s P dM dR PP INV MAL
dM=VM(:,1:qq)*diag(1./s(1:qq))*VD(:,1:qq)'*dR;
if INV.redu==1
    if INV.mitschicht==1
        dM=[P zeros(length(P(:,1)),length(PP(:,1)));PP]'*dM;
    else
        dM=P'*dM;
    end
end
nM=modelupdate(Model,dM,(INV.lolo==1));
figure(1);
draw3dmodel(nM,MAL);