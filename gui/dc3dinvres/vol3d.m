function varargout = vol3d(varargin)
% VOL3D Application M-file for vol3d.fig
%    FIG = VOL3D launch vol3d GUI.
%    VOL3D('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 24-Mar-2004 23:18:01

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    iconify(fig);
    
	if nargout > 0
		varargout{1} = fig;
	end
    global x y z M
    xm=x(1:end-1)+diff(x)/2;
    ym=y(1:end-1)+diff(y)/2;
    zm=z(1:end-1)+diff(z)/2;
    mi=log(min(M(:)));
    ma=log(max(M(:)));
    val=exp(mi+(ma-mi)*0.75);
    mal=1;
    while val>10, val=val/10;mal=mal*10; end
    val=rndig(val)*mal;
    set(handles.isoval,'String',num2str(val));
    set(handles.xslice,'String',sprintf('%.2f ',[min(xm) max(xm)]));
    set(handles.yslice,'String',sprintf('%.2f ',[min(ym) max(ym)]));
    set(handles.zslice,'String',sprintf('%.2f',max(zm)));
    draw3dgridmodel(handles);
    uiwait(fig);

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end

end

% --------------------------------------------------------------------
function draw3dgridmodel(handles)
global x y z M libmmfile MAL N
figure(5);
if ~isequal(libmmfile,4), set(5,'MenuBar','none','NumberTitle','off','Name','3D-Model'); end
iconify(5); % Testversion 
Xm=x(1:end-1)+diff(x)/2;
Ym=y(1:end-1)+diff(y)/2;
Zm=z(1:end-1)+diff(z)/2;
xstr=get(handles.xslice,'String');
if ~isempty(xstr), if(xstr=='.'), set(handles.xslice,'String',...
   sprintf('%.2f ',[Xm(1) Xm(end)])); end, end
ystr=get(handles.yslice,'String');
if ~isempty(ystr), if(ystr=='.'), set(handles.yslice,'String',...
   sprintf('%.2f ',[Ym(1) Ym(end)])); end, end
zstr=get(handles.zslice,'String');
if ~isempty(zstr), if(zstr=='.'), set(handles.zslice,'String',...
   sprintf('%.2f ',Zm(end))); end, end
if get(handles.smooth,'Value'),
    MM=smooth3(M,'gaussian');
else
    MM=M;
end
modelvolume(MM,x,y,z,str2num(get(handles.isoval,'String')),...
    str2num(get(handles.xslice,'String')),...
    str2num(get(handles.yslice,'String')),...
    str2num(get(handles.zslice,'String')));
asp=str2num(get(handles.high,'String'));
set(gca,'DataAspectRatio',[1 1 1/asp]);
view(str2num(get(handles.az,'String')),str2num(get(handles.el,'String')));
if isfield(MAL,'elec')&&any(MAL.elec)&isfield(N,'elec'),
    hold on;plot3(N.elec(:,1),N.elec(:,2),N.elec(:,3),'w.','MarkerSize',2);hold off
end
cb=colorbar;
set(cb,'DataAspectRatio',[64 1.5 1]);
yt=get(cb,'YTick');%ytl=get(cb,'YTickLabel');
ytl={};
% for i=1:length(yt), aa=num2str(round(10^yt(i)));ytl(i,1:length(aa))=aa; end
for i=1:length(yt), ytl{i}=num2str(rndig(10^yt(i))); end
set(cb,'YTickLabel',ytl);

% --------------------------------------------------------------------
function varargout = isoval_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = xslice_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = yslice_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = zslice_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = apply_Callback(h, eventdata, handles, varargin)
draw3dgridmodel(handles);

% --------------------------------------------------------------------
function varargout = close_Callback(h, eventdata, handles, varargin)
%uiresume(h);
uiresume(gcbf);
delete(gcbf);


% --- Executes during object creation, after setting all properties.
function az_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function az_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function el_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function el_Callback(hObject, eventdata, handles)

% --- Executes on button press in smooth.
function smooth_Callback(hObject, eventdata, handles)


% --- Executes on button press in export.
function export_Callback(hObject, eventdata, handles)
global datfile
if isempty(datfile),
    outfile='*.png';
else
    [path,name,ext]=fileparts(datfile); 
    outfile=strrep(datfile,ext,'-3d.png');
end
[outfile,iseps]=getgrfile(outfile); 
%[fname,pname]=uiputfile(outfile,'Save Figure as');
%if fname~=0,
%    outfile=fullfile(pname,fname);
%    [ff,pp,ee]=fileparts(outfile);
%    if strcmp(ee,''), outfile=[outfile '.png']; end
if ~isempty(outfile),
    message(['Exporting Figure to image ' outfile]);
    if iseps, % Testversion n. Zeile einklammern
        global libmmfile
        if ~isequal(libmmfile,4), errordlg('No eps export in Test version!');return; end
        print(5,'-depsc2',outfile);
    else exportpng(5,outfile); end
end 


% --- Executes during object creation, after setting all properties.
function high_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function high_Callback(hObject, eventdata, handles)


% --- Executes on button press in getazel.
function getazel_Callback(hObject, eventdata, handles)
if ishandle(5),
    set(figure(5),'MenuBar','none');
    [az,el]=view;
    set(handles.az,'String',num2str(az));    
    set(handles.el,'String',num2str(el));    
end
figure(handles.figure1);
