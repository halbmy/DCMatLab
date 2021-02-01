function message(text)

% MESSAGE - Outputs text into message box or on stdout
% message(text)

global output
if ishandle(output),
  ss=get(output,'String');
  ss{length(ss)+1}=text;
  unit=get(output,'Units');
  set(output,'String',ss,'Units','Characters');
  pos=get(output,'position');
  len=fix(pos(4))-2;
  if length(ss)>len, set(output,'ListboxTop',length(ss)-len); end
%   set(output,'Units',unit);
else
  fprintf('%s!!!\n',text);
end
    
