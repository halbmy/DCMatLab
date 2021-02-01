function pro2pdf(profile),

% PRO2PDF - Creates pdf-file displaying data in pro-file
% pro2pdf(profile)
% profile('test.pro') creates 'test.pro-pdf'

if nargin<1, error('Filename as input required!'); end
fid=fopen(profile,'r');
if fid<0, error('Error in opening file'); end
zeile=fgets(fid),
i=0;cmin=1e10;cmax=0;
while ischar(zeile),
    i=i+1;
    datei=sscanf(zeile,'%s %f %f %f %f',1);
    fprintf('%s\n',datei);
    switch check2dfile(datei),
        case 1,
            N=readres2dinvfile(datei);
        case 2,
            N=readinv2dfile(datei);
    end 
    cmin=min(cmin,min(N.r));
    cmax=max(cmax,max(N.r));
    NN{i}=N;
    name{i}=datei;
    zeile=fgets(fid);
end
fclose(fid);
texfile=[profile '.tex'];
fid=fopen(texfile,'w');
fprintf(fid,'\\documentclass[12pt,pdftex]{scrartcl}\n');
fprintf(fid,'\\usepackage{graphicx}\n');
fprintf(fid,'\\begin{document}\n');
fprintf(fid,'\\begin{center}\n');
% if cmin<1, cmin=1; end
% if cmax>1000, cmax=1000; end
mal.cauto=0;mal.cmin=cmin;mal.cmax=cmax;mal.log=(cmin>0);
showdata2d(NN{1},[],mal);
fprintf('resize figure window...\n');
pause
set(gca,'Units','pixels');
set(gcf,'Units','pixels');
for i=1:length(NN),
    fprintf('%d ',i);
    if mod(i,10)==0, fprintf('\n'); end
    showdata2d(NN{i},[],mal);
    set(gca,'FontSize',12);
    drawnow;
    po=get(gca,'position');
    set(gcf,'PaperSize',po(3:4),'PaperPositionMode','auto'); 
    if 0,
        ifile=[num2str(i) '.eps'];
        print(gcf,'-depsc2',ifile);
        dos(['epstopdf ' ifile])
        delete(ifile);
    else
        ifile=[num2str(i) '.png'];
        print(gcf,'-dpng','-r200',ifile);
    end
    fprintf(fid,'\\par %s\\\\\n',name{i});
    fprintf(fid,'\\includegraphics[height=0.3\\textheight]{%d}\n',i);
end
fprintf(fid,'\\end{center}\n');
fprintf(fid,'\\end{document}\n');
fclose(fid);
dos(['pdflatex -interaction=nonstopmode ' texfile]);