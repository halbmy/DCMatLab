function dc3drun(varargin)
global Model N S MAL INV FOR
load setup
if ~isequal(checklic3,4),
    fprintf('No valid license file found!\nExiting\n');
    return;
end
Err=struct('proz',3,'umin',100e-6,'curr',0.1);
zmax=0;nlay=11;xtra=1;INV.auto=2;Ms=[];outfile='';
if nargin==0, helptext;return; end
for i=1:length(varargin),
    argi=varargin{i};
    fprintf('Arg %d: ',i);
    if isequal(argi(1),'-'),
        ch=argi(2);rest=argi(3:end);
        switch argi(2),
            case 'h', helptext;
            case 'l',
                INV.lam=str2num(rest);
                fprintf('set lambda to %f\n',INV.lam);
            case 'O', INV.auto=1;INV.lam=1;
                fprintf('Set lambda optimization\n');
            case 's', [Ms.M,Ms.x,Ms.y,Ms.z]=modelimport3d(rest);
            case 'o', outfile=rest;
		fprintf('set outfile to %s\n',outfile);
            case 'z', zmax=str2num(rest);
                fprintf('Setting maximum z to %.1f\n',zmax);
            case 'n', nlay=num2str(rest);fprintf('using %d layers\n',nlay);
        end
    else
        datfile=argi;
        fprintf('processing data file %s\n',datfile);
        switch check3dfile(datfile),
            case 1,
                N=readres3dinvfile(datfile);
            case 2,
                N=readinv3dfile(datfile);
            otherwise
                error('File type unknown!');
        end 
        if ~isfield(N,'i'), N.i=Err.curr; end
        N.err=estimateerror(N,Err.proz,Err.umin,N.i); % 3% plus 100myV/100mA
        fprintf('Error estimate: %.1f%% + %dmyV/%dmA: min=%.1f%% max=%.1f%%\n',...
            Err.proz,round(Err.umin*1e6),round(Err.curr*1e3),min(N.err)*100,max(N.err)*100);
        if ~isempty(Ms), M=Ms;x=xs;z=zs;
        else
            Model.z=slimzparam(N,nlay,zmax);
            del=diff(unique(N.elec(:,1)));dx=min(del(find(del>0.1)));
            del=diff(unique(N.elec(:,2)));dy=min(del(find(del>0.1)));
            Model.x=min(N.elec(:,1))-xtra*dx:dx:max(N.elec(:,1))+xtra*dx;
            Model.y=min(N.elec(:,2))-xtra*dy:dy:max(N.elec(:,2))+xtra*dy;
            rhostart=median(N.r);Model.Bg=rhostart;
            Model.M=ones(length(Model.x)-1,length(Model.y)-1,length(Model.z)-1)*rhostart;
            fprintf('Model build of %dx%dx%d=%d cells, x=%.1f..(%.1f)..%.1f y=%.1f..(%.1f)..%.1f z=0..%.1f(%d layers)\n',...
                size(Model.M,1),size(Model.M,2),size(Model.M,3),prod(size(Model.M)),...
                min(Model.x),dx,max(Model.x),min(Model.y),dy,max(Model.y),max(Model.z),length(Model.z)-1);
        end
        S=getsens3d(datfile,Model,N);Cov=Model.M;Cov(:)=sum(abs(S));
        R=ones(size(N.r))*rhostart;
        M0=Model.M;
        RMSerr=rms(N.r,R,INV.lolo);
        CHIQ=chi2(N.r,R,N.err);
        fprintf('RMS = %.1f Chi^2 = %.1f\n',RMSerr,CHIQ);
        Mref=Model.M;
        running=1;
        while running,
            dR=log(N.r(:)-INV.lbound)-log(R(:)-INV.lbound);
            dM0=log(Model.M(:)-INV.lbound)-log(Mref(:)-INV.lbound);
            [dM,INV.lam]=slimminvers(S,dR,INV,N.err,Model,dM0);
            if (INV.auto==1), INV.auto=2; end
            newModel=modelupdate(Model,dM,2-INV.lolo,INV.lbound);
            oldR=R;
            R=mfdfwd3d(newModel,N,FOR);
            if isfield(INV,'linesearch')&&(INV.linesearch>0),
                tau=linesearch(N,oldR,R);
                if tau==0, message('tau=0, stopping inversion');break; end
                if tau>0.9, tau=1; end
                if tau<1, 
                    fprintf('line search parameter %.2f\n',tau);
                    newModel=modelupdate(Model,dM,2-INV.lolo,INV.lbound);
                    R=mfdfwd3d(newModel,N,FOR);
                end
            end
            RMSerr=[RMSerr rms(N.r,R,INV.lolo)];
            CHIQ=[CHIQ chi2(N.r,R,N.err)];
            fprintf('RMS = %.1f Chi^2 = %.1f\n',RMSerr(end),CHIQ(end));
            dchiq=CHIQ(end)-CHIQ(end-1);
            if isfield(INV,'lolo')&&(INV.lolo==0),
                ddr=(R-oldR-S*dM(:))/(dM(:)'*dM(:));ddr=ddr(:);
            else
                ddr=(log(R)-log(oldR)-S*dM(:))/(dM(:)'*dM(:));ddr=ddr(:);
            end
            for i=1:size(S,2), S(:,i)=S(:,i)+ddr*dM(i); end 
            message('Doing Broyden Update of Sensitivity');
            if (dchiq<0)|(length(CHIQ)<3),
                Model=newModel;
                running=(abs(dchiq)/CHIQ(end)>0.05);
                if running==0, message('Stopping (dCHI^2<5 percent)'); end
            else
                message('Going back to old model and stopping');
                running=0;
            end
            if CHIQ(end)<0.9, running=0; end
        end
        if isempty(outfile),
            [path,name,ext]=fileparts(datfile); 
            %outfile=strrep(datfile,ext,['-' num2str(INV.lam) '.mod']); 
            outfile=strrep(datfile,ext,'.mod'); 
        end
        modelexport3d(outfile,Model.M,Model.x,Model.y,Model.z);
        fprintf('exporting model to file %s\n',outfile);
        outfile='';
    end
end

function helptext
fprintf('DC3DRUN - Script version of DC3DInvRes\n');
fprintf('--------------------------------------\n');
fprintf('dc3drun <options> datafile1 <datafile2 ...>\n');
fprintf('Options: (switch-value-pairs are without blank!)\n');
fprintf('-h\t\t\tdisplay this help text\n');
fprintf('-llambda\tset regularization parameter\n');
fprintf('-O\t\t\toptimize regularization parameter\n');
fprintf('-ooutfile\tset name of model outputfile\n');
fprintf('-smodelfile\tset starting model\n');
fprintf('\n');
fprintf('\n');
