function isok=checklic2()

isok=0;if isunix, isok=4;return; end
file='dc2dinvres.lic';
fid=fopen(file,'r');
if fid>0,
%     s=fscanf(fid,'%s\r\n');
    s1=fgetl(fid);
    s1=fgetl(fid);
    s2=fgetl(fid);
    s3=fgetl(fid);
    fclose(fid);
    d1=double(s1)-48;
    d2=double(s2)-48;
    d3=double(s3)-48;
    eins=d2([1 2 3 5 8 13 21 34 35 36 6 7]);
    zwei=d2(22:33)/2;
    drei=double(convadd(getmacaddress))-48;
    ind=[1:2:35 36:-2:2];
    if isequal(d1,51-d2), isok=isok+1; end
    if isequal(d2(ind)+2,d3), isok=isok+1; end
    if isequal(eins,zwei), isok=isok+1; end
    if isequal(eins,drei), isok=isok+1; end
end    
