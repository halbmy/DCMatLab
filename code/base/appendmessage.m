function appendmessage(text)

% APPENDMESSAGE - Outputs text into message box or on stdout
% appendmessage(text)

global output
if ishandle(output),
  ss=get(output,'String');
  ss{end}=[ss{end} text];
  set(output,'String',ss,'Units','Characters');
  pos=get(output,'position');
  len=fix(pos(4))-2;
  if length(ss)>len, set(output,'ListboxTop',length(ss)-len); end
else
  fprintf('%s!!!\n',text);
end
    
