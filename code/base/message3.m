function message(text)

global output
if isempty(output)||(~ishandle(output)),
    fprintf('%s!!!\n',text);
else
    ss=get(output,'String');
    ss{length(ss)+1}=text;
    set(output,'String',ss);
    len=6;
    if length(ss)>len, set(output,'ListboxTop',length(ss)-len); end
end
%pause(0.1);
