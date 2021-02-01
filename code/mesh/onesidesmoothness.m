function C=onesidesmoothness(Mesh,po,usew)

% ONESIDESMOOTHNESS - Unstructured mesh smoothness operator (unsymmetric)
% C=onesidesmoothness(Mesh)
% C=onesidesmoothness(Mesh,zpower) % vertical smoothness using zpower
% C=onesidesmoothness(Mesh,zpower,1) % use zweight instead of zpower

if nargin<3, usew=0; end
if nargin<2, po=0; end
if Mesh.dim==2,
    blr=[Mesh.boundleft(:) Mesh.boundright(:)];
    bn=Mesh.bound;
    [i,j]=find(blr==0);
    blr(i,:)=[];bn(i,:)=[];
    nblr=size(blr,1);
    C=spalloc(nblr,Mesh.ncells,nblr*2);
    for i=1:nblr,
        if usew>0, 
            nv=diff(Mesh.node(bn(i,:),:));
            nv=nv/norm(nv);
%             val=po^(abs(nv(1))/usew);
            val=(1+po-abs(nv(1)))^usew;
        else
            %     ve=diff(Mesh.node(Mesh.bound(blr(i,:),:),:));
            ve=diff(Mesh.node(bn(i,:),:));
            val=(1-abs(ve(1))/norm(ve))^po; 
        end
        C(i,blr(i,:))=[-1 1]*val;    
    end
end
if Mesh.dim==3,
    C=spalloc(Mesh.ncells*4,Mesh.ncells,Mesh.ncells*8);
    for i=1:Mesh.ncells,        
        [i1,j1]=find(Mesh.cell==Mesh.cell(i,1));
        [i2,j2]=find(Mesh.cell==Mesh.cell(i,2));
        [i3,j3]=find(Mesh.cell==Mesh.cell(i,3));
        [i4,j4]=find(Mesh.cell==Mesh.cell(i,4));        
        ii=intersect(intersect(i1,i2),i3);
        if length(ii)>1, C((i-1)*4+1,ii)=[-1 1]; end
        ii=intersect(intersect(i1,i2),i4);
        if length(ii)>1, C((i-1)*4+2,ii)=[-1 1]; end
        ii=intersect(intersect(i1,i3),i4);
        if length(ii)>1, C((i-1)*4+3,ii)=[-1 1]; end
        ii=intersect(intersect(i2,i3),i4);
        if length(ii)>1, C((i-1)*4+4,ii)=[-1 1]; end        
    end
    C(sum(abs(C),2)==0,:)=[];
end