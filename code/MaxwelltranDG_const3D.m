function MaxwelltranDG_const3D%Total 2nd edge element for constant case: N=16
load('totalquadedgeelement3D.mat')
eta1=100;%eta=eta1/100
node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]/2+0.5;%/2;
elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7];
% elemH = label3(nodeH,elemH);        % label the mesh,  only used in adative refinement
for k=1:3
    [node,elem] = uniformbisect3(node,elem);%showmesh3(node,elem,[130,28],0.35);   
end
% n=5;clear node elem%Cube
% h = 2/n; [x y z] = meshgrid(-1:h:1,-1:h:1,-1:h:1);
% [cx cy cz] = meshgrid(-1+h/2:h:1-h/2,-1+h/2:h:1-h/2,-1+h/2:h:1-h/2);
% node(:,1) = [x(:); cx(:)]/2+.5*1;%[-0.5 0.5] 和[0 1] 计算时volume 稍不同
% node(:,2) = [y(:); cy(:)]/2+.5*1;
% node(:,3) = [z(:); cz(:)]/2+.5*1; elem = delaunayn(node);
% clear node elem%thickL
%  h = 1/1; [x y z] = meshgrid(-1:h:1,-1:h:1,-1:h:1);
% [cx cy cz] = meshgrid(-1+h/2:h:1-h/2,-1+h/2:h:1-h/2,-1+h/2:h:1-h/2);
% node(:,1) = [x(:); cx(:)]; node(:,2) = [y(:); cy(:)]; node(:,3) = [z(:); cz(:)];  
% nodedrop= node(:,1)<0&node(:,2)<0|node(:,3)<0; 
% node=node(~nodedrop,:);%thick L
% elem = delaunayn(node); 
% nodeHx=node(:,1); nodeHy=node(:,2);nodeHz=node(:,3); nodeHx=nodeHx(elem); nodeHy=nodeHy(elem);nodeHz=nodeHz(elem);
% elemdrop=all(nodeHx<=0&nodeHy<=0|nodeHz<=0,2);elem=elem(~elemdrop,:);  tetramesh(elem,node)
% elem = label3(node,elem);

bdFlag = setboundary3(node,elem,'Dirichlet');   % boundary faces
% clear    % L-shaped domain
% nodeH = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
% elemH = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7];elemH = label3(nodeH,elemH);        
% for k=1:1
%     [nodeH,elemH] = uniformbisect3(nodeH,elemH);%showmesh3(node,elem,[130,28],0.35);   
% end
%     nodeHx=nodeH(:,1); nodeHy=nodeH(:,2);nodeHz=nodeH(:,3); nodeHx=nodeHx(elemH); nodeHy=nodeHy(elemH);nodeHz=nodeHz(elemH);
%     elemdrop=all(nodeHx<=0&nodeHy<=0|nodeHz<=0,2);elemH=elemH(~elemdrop,:);
%     nodedrop= nodeH(:,1)<0&nodeH(:,2)<0|nodeH(:,3)<0; node=nodeH(~nodedrop,:);
%      [~,~,position]=unique(elemH(:)); nodenumloc1=1:size(nodeH,1);
%     elemH=nodenumloc1(position); elem=reshape(elemH,length(elemH)/4,4); clear nodenumloc1 position
% showmesh3(node,elem,[130,28],0.35);
% findnode3(nodeH);
% thick_L_mesh_local_refinement;%局部加密有厚度L区域
% for k=1:1
% %     [node,elem,~,bdFlag] = uniformbisect3(node,elem,[],bdFlag);%showmesh3(node,elem,[130,28],0.35);   
%     [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
% end

[elem,bdFlag] = sortelem3(elem,bdFlag);
%not suited for the longest-edge rule
%% Construct Data Structure
% [Dlambda,volume] = gradbasis3(node,elem);
 T= auxstructure3(elem(:,1:4)); NT=size(elem,1);%elem2face=T.elem2face;
 face=T.face;face2elem=T.face2elem; %Du = [dudx, dudy]; size(Du,1)=NT;行数
% [dudxc1,dudyc2,jumpdudxc1,jumpdudxc2,jumpdudyc1,jumpdudyc2,volume,Dlambda] = gradientconstant(node,elem,u,neighbor); 
v12 = node(face(:,2),:)-node(face(:,1),:);
v13 = node(face(:,3),:)-node(face(:,1),:);
normal = mycross(v12,v13,2);
facel = sqrt(sum(normal.^2,2))/2;

facet = uint32([elem(:,[2 4 3]);elem(:,[1 3 4]);elem(:,[1 4 2]);elem(:,[1 2 3])]);
v12 = node(facet(:,2),:)-node(facet(:,1),:);
v13 = node(facet(:,3),:)-node(facet(:,1),:);
normal = mycross(v12,v13,2);
v14 = node(elem(:,4),:)-node(elem(:,1),:);
volume = dot(normal(3*NT+1:4*NT,:),v14,2)/6;
Dlambda = zeros(NT,3,4);
Dlambda(1:NT,:,1) = normal(1:NT,:)./[6*volume, 6*volume, 6*volume];
Dlambda(1:NT,:,2) = normal(NT+1:2*NT,:)./[6*volume, 6*volume, 6*volume];
Dlambda(1:NT,:,3) = normal(2*NT+1:3*NT,:)./[6*volume, 6*volume, 6*volume];
Dlambda(1:NT,:,4) = normal(3*NT+1:4*NT,:)./[6*volume, 6*volume, 6*volume];
idx = (volume<0); 
volume(idx,:) = -volume(idx,:);
for i=1:4;normal(NT*(i-1)+1:NT*i,:)=-normal(NT*(i-1)+1:NT*i,:)./repmat(sqrt(sum(normal(NT*(i-1)+1:NT*i,:).^2,2)).*(-1).^idx,1,3);end
% dir=accumarray([face2elem(:,1),face2elem(:,3)],1,[NT,4]);%取边的第一个单元的外法向基函数（囊括了边界边）
% inedgejudge=face2elem(:,1)~=face2elem(:,2);%内部边逻辑向量
% dir=dir+accumarray([face2elem(inedgejudge,2),face2elem(inedgejudge,4)],-1,[NT,4]);
% td=[ones(NT,18),repmat(dir,1,3)]; d=zeros(NT,30,30);
% for i=1:30;
%     for j=1:30
%         d(:,i,j)=td(:,i).*td(:,j);
%     end
% end
B=sym('[B11 B12 B13;B21 B22 B23;B31 B32 B33]');
%det(B)=-volume*(-1).^idx*6, %volume*(-1).^idx是用12 13 14 vector计算的行列式/6,
%det(B):41 42 43  vector计算的行列式, 两者正好相差-1
Dl=sym('[D1l1 D2l1 D3l1;D1l2 D2l2 D3l2;D1l3 D2l3 D3l3]');
u=Dl.'*basis;%u(1:3,1:30)
c1u=-B*curlbasis;%c1u(1:3,1:30), not divide by -|B|=-|A41 A42 A43|=|A12 A13 A14|
c2u=sym(zeros(3,30)); 
for i=1:30;
    Dcurlbasis=[diff(curlbasis(:,i),'lam1') diff(curlbasis(:,i),'lam2') diff(curlbasis(:,i),'lam3')];
    Dxcurlu=-B*Dcurlbasis*Dl;%not divide by -|B|
    c2u(:,i)=[Dxcurlu(3,2)-Dxcurlu(2,3);Dxcurlu(1,3)-Dxcurlu(3,1);Dxcurlu(2,1)-Dxcurlu(1,2)];
end
B1=node(elem(:,1),:)-node(elem(:,4),:);  B2=node(elem(:,2),:)-node(elem(:,4),:); B3=node(elem(:,3),:)-node(elem(:,4),:); 
% B11=B1(:,1); B21=B1(:,2); B31=B1(:,3); B12=B2(:,1); B22=B2(:,2); B32=B2(:,3);
% B13=B3(:,1); B23=B3(:,2); B33=B3(:,3); 
Dl1=Dlambda(:,:,1); Dl2=Dlambda(:,:,2); Dl3=Dlambda(:,:,3);
% D1l1=Dl1(:,1); D2l1=Dl1(:,2); D3l1=Dl1(:,3); 
% D1l2=Dl2(:,1); D2l2=Dl2(:,2); D3l2=Dl2(:,3); D1l3=Dl3(:,1); D2l3=Dl3(:,2); D3l3=Dl3(:,3);
% locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
% locFace = [2 3 4; 1 3 4; 1 2 4; 1 2 3];
% locBasesIdx = [1 2 0; 1 3 0; 1 4 0; 2 3 0; 2 4 0; 3 4 0; ... % phi
%                1 2 0; 1 3 0; 1 4 0; 2 3 0; 2 4 0; 3 4 0; ... % psi
%                3 2 4; 3 1 4; 2 1 4; 2 1 3; ...
%                4 2 3; 4 1 3; 4 1 2; 3 1 2]; % chi
% See <matlab:ifem('Maxwell2doc') Maxwell2doc> for details.            
% construct elem2edge and elem2face 
[elem2edge,edge,elem2edgeSign] = dof3edge(elem);
[elem2face,face] = dof3face(elem); 
NT = size(elem,1);  %N = size(node,1);
NE = size(edge,1); NF = size(face,1); Ndof = 3*(NE + NF);
elem2dof = [elem2edge elem2edge+NE elem2edge+2*NE elem2face+3*NE elem2face+3*NE+NF elem2face+3*NE+2*NF];
face2edge = zeros(NF,3,'int32');
% face is given by auxstructure3 and is sorted according to global indices.
% locFace = [2 3 4; 1 3 4; 1 2 4; 1 2 3];
% face2edge is used to compute uI. So it is consistent with the local index
% system in edgeinterpolate2, i.e., if face is (i,j,k) with i<j<k, then the
% three edges are [i j], [i k], [j k]. 
face2edge(elem2face(:,1),:) = elem2edge(:,[4 5 6]);%面的边按升序排,此程序没用到face
face2edge(elem2face(:,2),:) = elem2edge(:,[2 3 6]);
face2edge(elem2face(:,3),:) = elem2edge(:,[1 3 5]);
face2edge(elem2face(:,4),:) = elem2edge(:,[1 2 4]);

[lambda,weight] = quadpts3(4);lambda1=lambda(:,1)';lambda2=lambda(:,2)';lambda3=lambda(:,3)'; lambda4=lambda(:,4)'; 
lam1=repmat(lambda1,NT,1); lam2=repmat(lambda2,NT,1); lam3=repmat(lambda3,NT,1);
weight=weight(:)';
nQuadT=size(weight,2); %x0=(node(elem(:,1),1)+node(elem(:,2),1)+node(elem(:,3),1)+node(elem(:,4),1))/4; 
% y0=(node(elem(:,1),2)+node(elem(:,2),2)+node(elem(:,3),2)+node(elem(:,4),2))/4;
% z0=(node(elem(:,1),3)+node(elem(:,2),3)+node(elem(:,3),3)+node(elem(:,4),3))/4;
% inv([16 x y;x 16 z;y z 14]-eye(3))
x=kron(node(elem(:,1),1),lambda1)+kron(node(elem(:,2),1),lambda2)+kron(node(elem(:,3),1),lambda3)+kron(node(elem(:,4),1),lambda4);
y=kron(node(elem(:,1),2),lambda1)+kron(node(elem(:,2),2),lambda2)+kron(node(elem(:,3),2),lambda3)+kron(node(elem(:,4),2),lambda4);
z=kron(node(elem(:,1),3),lambda1)+kron(node(elem(:,2),3),lambda2)+kron(node(elem(:,3),3),lambda3)+kron(node(elem(:,4),3),lambda4);
% nv11=(z^2 - 195)/(13*x^2 - 2*x*y*z + 15*y^2 + 15*z^2 - 2925);
% nv12=(13*x - y*z)/(13*x^2 - 2*x*y*z + 15*y^2 + 15*z^2 - 2925);
% nv13=(15*y - x*z)/(13*x^2 - 2*x*y*z + 15*y^2 + 15*z^2 - 2925);
% nv21=nv12; nv22=(y^2 - 195)/(13*x^2 - 2*x*y*z + 15*y^2 + 15*z^2 - 2925);
% nv23=(15*z - x*y)/(13*x^2 - 2*x*y*z + 15*y^2 + 15*z^2 - 2925);
% nv31=nv13; nv32=nv23; nv33=(x^2 - 225)/(13*x^2 - 2*x*y*z + 15*y^2 + 15*z^2 - 2925);
nd=16;%n=8+xt-yt;
% y0=repmat(y0,1,nQuadT); x0=repmat(x0,1,nQuadT); z0=repmat(z0,1,nQuadT);
% xt=xt-x0; yt=yt-y0; zt=zt-z0;
B11=repmat(B1(:,1),1,nQuadT); B21=repmat(B1(:,2),1,nQuadT); B31=repmat(B1(:,3),1,nQuadT); B12=repmat(B2(:,1),1,nQuadT); 
B22=repmat(B2(:,2),1,nQuadT); B32=repmat(B2(:,3),1,nQuadT);
B13=repmat(B3(:,1),1,nQuadT); B23=repmat(B3(:,2),1,nQuadT); B33=repmat(B3(:,3),1,nQuadT); 
D1l1=repmat(Dl1(:,1),1,nQuadT); D2l1=repmat(Dl1(:,2),1,nQuadT); D3l1=repmat(Dl1(:,3),1,nQuadT); 
D1l2=repmat(Dl2(:,1),1,nQuadT); D2l2=repmat(Dl2(:,2),1,nQuadT); D3l2=repmat(Dl2(:,3),1,nQuadT); 
D1l3=repmat(Dl3(:,1),1,nQuadT); D2l3=repmat(Dl3(:,2),1,nQuadT); D3l3=repmat(Dl3(:,3),1,nQuadT);
% lambda=zeros(NT,nQuad,3);
weight=repmat(weight,NT,1); %size(wd)
phi=zeros(NT,nQuadT,30,3); curlphi=phi; curl2phi=phi;  %bt=zeros(NT,nQuadT,3);
for j=1:30
    phi(:,:,j,1)=eval(u(1,j)); phi(:,:,j,2)=eval(u(2,j));
    phi(:,:,j,3)=eval(u(3,j));
    curlphi(:,:,j,1)=eval(c1u(1,j)); curlphi(:,:,j,2)=eval(c1u(2,j));
    curlphi(:,:,j,3)=eval(c1u(3,j));
    curl2phi(:,:,j,1)=eval(c2u(1,j)); curl2phi(:,:,j,2)=eval(c2u(2,j));
    curl2phi(:,:,j,3)=eval(c2u(3,j));
end
%% Triangulation output
% T = struct('edge',edge,'face',face,'face2edge',face2edge, ...
%            'elem',elem,'bdFlag',bdFlag);
ii = zeros(900*NT,1); jj = ii; iit = zeros(465*NT,1); jjt = iit; 
index = 0; indext = 0; sA=jj; sAt=sA; sB=iit; sD=sB; sDt=sB; sAc=sB;
% [lambda,w] = quadpts3(3); % quadrature order is 3
% nQuad = size(lambda,1);
for i = 1:30
    for j = 1:30
        ii(index+1:index+NT) = double(elem2dof(:,i));
        jj(index+1:index+NT) = double(elem2dof(:,j));
        %         i1 = locBasesIdx(i,1); i2 = locBasesIdx(i,2); i3 = locBasesIdx(i,3);
        t=(curl2phi(:,:,i,1).*phi(:,:,j,1)+curl2phi(:,:,i,2).*phi(:,:,j,2)+curl2phi(:,:,i,3).*phi(:,:,j,3));
        Aij=sum(t.*weight,2)/6.*(-1).^idx; Atij=sum((nd-1).^(-1).*t.*weight,2)/6.*(-1).^idx;%.*d(:,i,j)
        sA(index+1:index+NT) = Aij; sAt(index+1:index+NT) = Atij; index = index + NT;
        if j>=i
            iit(indext+1:indext+NT) = double(elem2dof(:,i));
            jjt(indext+1:indext+NT) = double(elem2dof(:,j));
            Bij=sum((curlphi(:,:,i,1).*curlphi(:,:,j,1)+curlphi(:,:,i,2).*curlphi(:,:,j,2)+curlphi(:,:,i,3).*curlphi(:,:,j,3)).*weight,2)/36./volume;
            %          A=A+sparse(elem2dof(:,i),elem2dof(:,j), sum(curlphi(:,:,i).*curlphi(:,:,j).*weight,2).*volume,Ndof,Ndof);
            t=(phi(:,:,i,1).*phi(:,:,j,1)+phi(:,:,i,2).*phi(:,:,j,2)+phi(:,:,i,3).*phi(:,:,j,3));
            Dij=sum(t.*weight,2).*volume; Dtij=sum((nd-1).^(-1).*t.*weight,2).*volume;
            Acij=sum((nd-1).^(-1).*(curl2phi(:,:,i,1).*curl2phi(:,:,j,1)+curl2phi(:,:,i,2).*curl2phi(:,:,j,2)+curl2phi(:,:,i,3).*curl2phi(:,:,j,3)).*weight,2)/36./volume;
            sB(indext+1:indext+NT) = Bij; sD(indext+1:indext+NT) = Dij;
            sDt(indext+1:indext+NT) = Dtij;
            sAc(indext+1:indext+NT) = Acij;
            indext = indext + NT;
        end
        %          M=M+sparse(elem2dof(:,i),elem2dof(:,j), sum((phi(:,:,i,1).*phi(:,:,j,1)+phi(:,:,i,2).*phi(:,:,j,2)).*weight,2).*volume,Ndof,Ndof);
    end
end
clear curlBasis_i curlBasis_j basis_i basis_j
A = sparse(ii,jj,sA,Ndof,Ndof); At = sparse(jj,ii,sAt,Ndof,Ndof); 
diagIdx = (iit == jjt);   upperIdx = ~diagIdx;
iid=iit(diagIdx); iiu=iit(upperIdx); jjd=jjt(diagIdx); jju=jjt(upperIdx);
D = sparse(iid,jjd,sD(diagIdx),Ndof,Ndof); Dt = sparse(jjd,iid,sDt(diagIdx),Ndof,Ndof);
DU = sparse(iiu,jju,sD(upperIdx),Ndof,Ndof);DtU = sparse(iiu,jju,sDt(upperIdx),Ndof,Ndof);
B = sparse(iid,jjd,sB(diagIdx),Ndof,Ndof); BU = sparse(iiu,jju,sB(upperIdx),Ndof,Ndof);% AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Ndof,Ndof);
B = B + BU + BU'; D = D + DU + DU'; Dt = Dt + DtU + DtU';
Ac = sparse(iid,jjd,sAc(diagIdx),Ndof,Ndof); AcU = sparse(iiu,jju,sAc(upperIdx),Ndof,Ndof);
Ac = Ac + AcU + AcU';
% M = sparse(ii(diagIdx),jj(diagIdx),sM(diagIdx),Ndof,Ndof);
% MU = sparse(ii(upperIdx),jj(upperIdx),sM(upperIdx),Ndof,Ndof);
% M = M + MU + MU';
% clear AU MU
 
[lambda,w] = quadpts(3);lambda1=lambda(:,1)';lambda2=lambda(:,2)';lambda3=lambda(:,3)'; 
lambda1=repmat(lambda1,NF,1); lambda2=repmat(lambda2,NF,1); lambda3=repmat(lambda3,NF,1);
lam1=zeros(size(lambda1)); lam2=lam1; lam3=lam1;
w=w(:)';nQuad=size(w,2);
    w=repmat(w,NF,1); %x0=x0';%x=repmat(x',NT,1);
% x(:,:,1)=kron(node(elem(:,2),1),x0)+kron(node(elem(:,3),1),1-x0);
% y(:,:,1)=kron(node(elem(:,2),2),x0)+kron(node(elem(:,3),2),1-x0);
% x(:,:,2)=kron(node(elem(:,3),1),x0)+kron(node(elem(:,1),1),1-x0);
% y(:,:,2)=kron(node(elem(:,3),2),x0)+kron(node(elem(:,1),2),1-x0);
% x(:,:,3)=kron(node(elem(:,1),1),x0)+kron(node(elem(:,2),1),1-x0);
% y(:,:,3)=kron(node(elem(:,1),2),x0)+kron(node(elem(:,2),2),1-x0);
% nx(:,:,1)=repmat(normal(:,1,1),1,nQuad); ny(:,:,1)=repmat(normal(:,2,1),1,nQuad);%第一条边的沿x y方向的单位法向
% nx(:,:,2)=repmat(normal(:,1,2),1,nQuad); ny(:,:,2)=repmat(normal(:,2,2),1,nQuad);%
% nx(:,:,3)=repmat(normal(:,1,3),1,nQuad); ny(:,:,3)=repmat(normal(:,2,3),1,nQuad);
% xx=kron(node(face(:,1),1),lambda1)+kron(node(face(:,2),1),lambda2)+kron(node(face(:,3),1),lambda3);
% yy=kron(node(face(:,1),2),lambda1)+kron(node(face(:,2),2),lambda2)+kron(node(face(:,3),2),lambda3);
% zz=kron(node(face(:,1),3),lambda1)+kron(node(face(:,2),3),lambda2)+kron(node(face(:,3),3),lambda3);
clear phi;phi{1}=zeros(NF,nQuad,30,3); phi{2}=phi{1}; phi{3}=phi{1}; curlphi=phi;curl2phi=phi;
for k=1:2
    %     xtt=(node(elem(face2elem(:,k),1),1)+node(elem(face2elem(:,k),2),1)+node(elem(face2elem(:,k),3),1)+node(elem(face2elem(:,k),4),1))/4;
    %     ytt=(node(elem(face2elem(:,k),1),2)+node(elem(face2elem(:,k),2),2)+node(elem(face2elem(:,k),3),2)+node(elem(face2elem(:,k),4),2))/4;
    %     ztt=(node(elem(face2elem(:,k),1),3)+node(elem(face2elem(:,k),2),3)+node(elem(face2elem(:,k),3),3)+node(elem(face2elem(:,k),4),3))/4;
    t=face2elem(:,k);
    B1=node(elem(t,1),:)-node(elem(t,4),:);  B2=node(elem(t,2),:)-node(elem(t,4),:); B3=node(elem(t,3),:)-node(elem(t,4),:);
    B11=repmat(B1(:,1),1,nQuad); B21=repmat(B1(:,2),1,nQuad); B31=repmat(B1(:,3),1,nQuad); B12=repmat(B2(:,1),1,nQuad);
    B22=repmat(B2(:,2),1,nQuad); B32=repmat(B2(:,3),1,nQuad);
    B13=repmat(B3(:,1),1,nQuad); B23=repmat(B3(:,2),1,nQuad); B33=repmat(B3(:,3),1,nQuad);
    Dl1=Dlambda(t,:,1); Dl2=Dlambda(t,:,2); Dl3=Dlambda(t,:,3);
    D1l1=repmat(Dl1(:,1),1,nQuad); D2l1=repmat(Dl1(:,2),1,nQuad); D3l1=repmat(Dl1(:,3),1,nQuad);
    D1l2=repmat(Dl2(:,1),1,nQuad); D2l2=repmat(Dl2(:,2),1,nQuad); D3l2=repmat(Dl2(:,3),1,nQuad);
    D1l3=repmat(Dl3(:,1),1,nQuad); D2l3=repmat(Dl3(:,2),1,nQuad); D3l3=repmat(Dl3(:,3),1,nQuad);
    for i=1:4
        t=face2elem(:,k+2)==i;
        switch i
            case 1
                lam1(t,:)=0; lam2(t,:)=lambda1(t,:); lam3(t,:)=lambda2(t,:);
            case 2
                lam2(t,:)=0; lam1(t,:)=lambda1(t,:); lam3(t,:)=lambda2(t,:);
            case 3
                lam3(t,:)=0; lam1(t,:)=lambda1(t,:); lam2(t,:)=lambda2(t,:);
            case 4
                lam1(t,:)=lambda1(t,:); lam2(t,:)=lambda2(t,:); lam3(t,:)=lambda3(t,:);
        end
    end
    % xtt=xx-repmat(xtt,1,nQuad); ytt=yy-repmat(ytt,1,nQuad); ztt=zz-repmat(ztt,1,nQuad);
    %     h=power(volume(face2elem(:,k)),1/3); h=repmat(h,1,nQuad);    x=xtt;y=ytt;z=ztt;
    vol=repmat(volume(face2elem(:,k)).*(-1).^idx(face2elem(:,k)),1,nQuad);
    for j=1:30
        phi{k}(:,:,j,1)=eval(u(1,j)); phi{k}(:,:,j,2)=eval(u(2,j));
        phi{k}(:,:,j,3)=eval(u(3,j));
        curlphi{k}(:,:,j,1)=eval(c1u(1,j))./vol/6; curlphi{k}(:,:,j,2)=eval(c1u(2,j))./vol/6;% divide by |B|
        curlphi{k}(:,:,j,3)=eval(c1u(3,j))./vol/6;
        curl2phi{k}(:,:,j,1)=eval(c2u(1,j))./vol/6; curl2phi{k}(:,:,j,2)=eval(c2u(2,j))./vol/6;
        curl2phi{k}(:,:,j,3)=eval(c2u(3,j))./vol/6;
    end
end
% direction=accumarray([face2elem(:,1),face2elem(:,3)],1,[NT,3]);%取边的第一个单元的外法向基函数（囊括了边界边）
% infacejudge=face2elem(:,1)~=face2elem(:,2);%内部边逻辑向量
% direction=direction+accumarray([face2elem(infacejudge,2),face2elem(infacejudge,4)],-1,[NT,3]);%除边界边外取边的第二个单元的内法向基函数
nx=zeros(NF,nQuad);ny=nx; nz=nx;
for i=1:4
    t=face2elem(face2elem(:,3)==i,1);
    nx(face2elem(:,3)==i,:)=repmat(normal(t+NT*(i-1),1),1,nQuad); ny(face2elem(:,3)==i,:)=repmat(normal(t+NT*(i-1),2),1,nQuad);
    nz(face2elem(:,3)==i,:)=repmat(normal(t+NT*(i-1),3),1,nQuad);
end
t=face2elem(:,1)==face2elem(:,2); %边界边
index=0; Nt=sum(t); Nnt=NF-Nt; indext=0;
sAnt=zeros(4*900*Nnt,1); sAt=zeros(900*Nt,1);  mm=sAnt;  nn=mm; sBnt=mm; mmt=sAt;  nnt=mmt; sBt=sAt;
sAAnt=sAnt; sAAt=sAt;
% sAt=zeros(9*dofl^2*Nt,1); mmt=sAt;  nnt=mmt;
for m=1:30%若组装时分成两个单元，组装时间还可以优化
    for n=1:30
%         A1=sum(0.5*curl3phi(:,:,m).*(phi(:,:,n,1).*ny-phi(:,:,n,2).*nx).*w,2).*facel;%单元自身基函数
%         A2=sum(0.5*curl3phi(:,:,n).*(phi(:,:,m,1).*ny-phi(:,:,m,2).*nx).*w,2).*facel;
        A1={0,0;0,0}; BB=A1;A2=A1;% B=A1;C=A1;%curl2uv={0,0};curl3uv={0,0};
%         Amn=sum(0.5*curl3phi(:,:,m).*(phi(:,:,n,1).*ny(:,:,i)-phi(:,:,n,2).*nx(:,:,i)).*w,2).*elem2facel(:,i);%相邻单元之间基函数（方法的优势）
%          Amn=Amn+sum(0.5*curl3phi(:,:,n).*(phi(:,:,m,1).*ny(:,:,i)-phi(:,:,m,2).*nx(:,:,i)).*w,2).*elem2facel(:,i);
        for k1=1:2%1-2号单元
            for k2=1:2
%         A1{k1,k2}=0.5*sum((curl2phi{k1}(:,:,m,1).*(curlphi{k2}(:,:,n,2).*nz-curlphi{k2}(:,:,n,3).*ny)...
%             +curl2phi{k1}(:,:,m,2).*(curlphi{k2}(:,:,n,3).*nx-curlphi{k2}(:,:,n,1).*nz)...
%             +curl2phi{k1}(:,:,m,3).*(curlphi{k2}(:,:,n,1).*ny-curlphi{k2}(:,:,n,2).*nx)+...
%            curl3phi{k1}(:,:,m,1).*(phi{k2}(:,:,n,2).*nz-phi{k2}(:,:,n,3).*ny)+...
%            curl3phi{k1}(:,:,m,2).*(phi{k2}(:,:,n,3).*nx-phi{k2}(:,:,n,1).*nz)+...
%            curl3phi{k1}(:,:,m,3).*(phi{k2}(:,:,n,1).*ny-phi{k2}(:,:,n,2).*nx))...
%             .*w,2).*facel;%单元自身基函数
        A1{k1,k2}=0.5*sum((phi{k1}(:,:,m,1).*(curlphi{k2}(:,:,n,2).*nz-curlphi{k2}(:,:,n,3).*ny)...
            +phi{k1}(:,:,m,2).*(curlphi{k2}(:,:,n,3).*nx-curlphi{k2}(:,:,n,1).*nz)...
            +phi{k1}(:,:,m,3).*(curlphi{k2}(:,:,n,1).*ny-curlphi{k2}(:,:,n,2).*nx))...
            .*w,2).*facel;%.*td(face2elem(:,k1),m).*td(face2elem(:,k2),n);
        A2{k1,k2}=0.5*sum((curl2phi{k1}(:,:,m,1).*(curlphi{k2}(:,:,n,2).*nz-curlphi{k2}(:,:,n,3).*ny)...
            +curl2phi{k1}(:,:,m,2).*(curlphi{k2}(:,:,n,3).*nx-curlphi{k2}(:,:,n,1).*nz)...
            +curl2phi{k1}(:,:,m,3).*(curlphi{k2}(:,:,n,1).*ny-curlphi{k2}(:,:,n,2).*nx))...
            .*w,2).*facel;
%         curl2uv{k2}=sum(curlphi{k2}(:,:,n).*(-curl2uexact(:,:,1).*ny+curl2uexact(:,:,2).*nx).*w,2).*facel;
%         curl3uv{k2}=sum(curl3uexact.*(phi{k2}(:,:,n,1).*ny(:,:,i)-phi{k2}(:,:,n,2).*nx(:,:,i)).*w,2).*facel;
%         A2{k1,k2}=0.5*sum((curl2phi{k2}(:,:,n,1).*(curlphi{k1}(:,:,m,2).*nz-curlphi{k1}(:,:,m,3).*ny)...
%             +curl2phi{k2}(:,:,n,2).*(curlphi{k1}(:,:,m,3).*nx-curlphi{k1}(:,:,m,1).*nz)...
%             +curl2phi{k2}(:,:,n,3).*(curlphi{k1}(:,:,m,1).*ny-curlphi{k1}(:,:,m,2).*nx)+...
%            curl3phi{k2}(:,:,n,1).*(phi{k1}(:,:,m,2).*nz-phi{k1}(:,:,m,3).*ny)+...
%            curl3phi{k2}(:,:,n,2).*(phi{k1}(:,:,m,3).*nx-phi{k1}(:,:,m,1).*nz)+...
%            curl3phi{k2}(:,:,n,3).*(phi{k1}(:,:,m,1).*ny-phi{k1}(:,:,m,2).*nx))...
%             .*w,2).*facel;%单元自身基函数
%         A1=A1+sum(curlphi(:,:,m).*(0.5*phi(:,:,n,1).*ny-0.5*phi(:,:,n,2).*nx).*w,2).*facel;%单元自身基函数
%         A2=A2+sum(curlphi(:,:,n).*(0.5*phi(:,:,m,1).*ny-0.5*phi(:,:,m,2).*nx).*w,2).*facel;
        BB{k1,k2}=sum(((curlphi{k1}(:,:,m,2).*nz-curlphi{k1}(:,:,m,3).*ny).*(curlphi{k2}(:,:,n,2).*nz-curlphi{k2}(:,:,n,3).*ny)...
            +(curlphi{k1}(:,:,m,3).*nx-curlphi{k1}(:,:,m,1).*nz).*(curlphi{k2}(:,:,n,3).*nx-curlphi{k2}(:,:,n,1).*nz)...
            +(curlphi{k1}(:,:,m,1).*ny-curlphi{k1}(:,:,m,2).*nx).*(curlphi{k2}(:,:,n,1).*ny-curlphi{k2}(:,:,n,2).*nx)).*w,2).*facel.^(1/2).*eta1;%单元自身基函数
%         B=zeros(size(A1));
%         C{k1,k2}=sum(((phi{k1}(:,:,m,2).*nz-phi{k1}(:,:,m,3).*ny).*(phi{k2}(:,:,n,2).*nz-phi{k2}(:,:,n,3).*ny)...
%             +(phi{k1}(:,:,m,3).*nx-phi{k1}(:,:,m,1).*nz).*(phi{k2}(:,:,n,3).*nx-phi{k2}(:,:,n,1).*nz)...
%             +(phi{k1}(:,:,m,1).*ny-phi{k1}(:,:,m,2).*nx).*(phi{k2}(:,:,n,1).*ny-phi{k2}(:,:,n,2).*nx)).*w,2).*facel.^(-1/2).*eta2*P^6;
%         C=sum((phi(:,:,m,1).*ny-phi(:,:,m,2).*nx).*(phi(:,:,n,1).*ny-phi(:,:,n,2).*nx).*w,2).*eta1;
            end
        end
%         sAnt(index+1:index+4*Nnt)=[A1{1,1}(~t)+A2{1,1}(~t)+B{1,1}(~t)+C{1,1}(~t); -A1{2,2}(~t)-A2{2,2}(~t)+B{2,2}(~t)+C{2,2}(~t);...
%              -A1{1,2}(~t)+A2{1,2}(~t)-B{1,2}(~t)-C{1,2}(~t); A1{2,1}(~t)-A2{2,1}(~t)-B{2,1}(~t)-C{2,1}(~t)]; 
        sAnt(index+1:index+4*Nnt)=[A1{1,1}(~t); -A1{2,2}(~t);-A1{1,2}(~t); A1{2,1}(~t)]; 
        sAAnt(index+1:index+4*Nnt)=[A2{1,1}(~t); -A2{2,2}(~t);-A2{1,2}(~t); A2{2,1}(~t)]; 
        sBnt(index+1:index+4*Nnt)=[BB{1,1}(~t); BB{2,2}(~t);-BB{1,2}(~t); -BB{2,1}(~t)]; 
        mm(index+1:index+4*Nnt) = elem2dof(face2elem(~t,[1 2 1 2]),m);
        nn(index+1:index+4*Nnt) = elem2dof(face2elem(~t,[1 2 2 1]),n);
        sAt(indext+1:indext+Nt)=2*A1{1,1}(t); sAAt(indext+1:indext+Nt)=2*A2{1,1}(t); sBt(indext+1:indext+Nt)=BB{1,1}(t);
        mmt(indext+1:indext+Nt) = elem2dof(face2elem(t,1),m);
        nnt(indext+1:indext+Nt) = elem2dof(face2elem(t,1),n);
        index=index+4*Nnt; indext=indext+Nt;
    end
end
AE =  sparse([mm;mmt],[nn;nnt],[sAnt;sAt],Ndof,Ndof); BE =  sparse([mm;mmt],[nn;nnt],[sBnt;sBt],Ndof,Ndof);
AAE =  sparse([mm;mmt],[nn;nnt],[sAAnt;sAAt],Ndof,Ndof);
%% Boundary conditions
%% Part 1: Find Dirichlet dof and modify the matrix
% Find Dirichlet boundary dof: fixedDof
    % Dirichlet boundary condition only
%     bdFlag = setboundary3(elem,1);
isBdDof = false(Ndof,1);

    %% Dirichlet boundary condition on edge dofs
    % Find boundary faces, edges and nodes
    isBdFace = false(NF,1);
    isBdFace(elem2face(bdFlag(:,1) == 1,1)) = true;
    isBdFace(elem2face(bdFlag(:,2) == 1,2)) = true;
    isBdFace(elem2face(bdFlag(:,3) == 1,3)) = true;
    isBdFace(elem2face(bdFlag(:,4) == 1,4)) = true;
%     boundaryFace = face(isBdFace,:);
    bdFace2edge = face2edge(isBdFace,:);
    isBdEdge = false(NE,1);
    isBdEdge(bdFace2edge(:)) = true;
    edgeBdDof = [find(isBdEdge); NE + find(isBdEdge); 2*NE + find(isBdEdge)];
%     bdEdge = edge(isBdEdge,:);
%     isBdNode(bdEdge) = true;
%     bdNode = find(isBdNode);
    faceBdDof = 3*NE + [find(isBdFace); NF+find(isBdFace); 2*NF+find(isBdFace)];
%     edgeIdxMap = zeros(NE,1);
%     edgeIdxMap(isBdEdge) = 1:size(bdEdge,1);
%     bdFace2edge = edgeIdxMap(bdFace2edge);
    isBdDof(edgeBdDof) = true;  
    isBdDof(faceBdDof) = true; freedof=[~isBdDof;~isBdDof];%freedof=[~isBdDof;true(Ndof,1);~isBdDof];
    tic;
    AA=[BE/100+Ac+(AAE'+AAE)/15,sparse(Ndof,Ndof);sparse(Ndof,Ndof), speye(Ndof)];
    BB=[At+At'+B+(AE'+AE)/15, -Dt-D; speye(Ndof),sparse(Ndof,Ndof)];
    [eigf,eigv]=eigs(AA(freedof,freedof),BB(freedof,freedof),7,5)%pause at this line
    diag(eigv.^.5)
    for i=1:length(eigv)
        indicator(i,:)=abs(  eigf(1:end/2,i)'*(-At(~isBdDof,~isBdDof)-At(~isBdDof,~isBdDof)'+...
            2*eigv(i,i)*(Dt(~isBdDof,~isBdDof)+D(~isBdDof,~isBdDof))-AE(~isBdDof,~isBdDof)/15-AE(~isBdDof,~isBdDof)'/15)*eigf(1:end/2,i)...
        /(eigf(1:end/2,i)'*B(~isBdDof,~isBdDof)*eigf(1:end/2,i))-1  );
    end
    indicator
    toc