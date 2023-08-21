function MaxwelltranDG_const2D%Total 2nd edge element for constant case: index of refraction n=16
load('totalquadedgeelement.mat')
eta1=100;%eta=eta1/100
node=[-3^0.5/2,-0.5; 3^0.5/2,-0.5; 0 1]; elem=[1 2 3]; %triangle
 node = [-1 1;-1 0;0 1;0 0;1 1;1 0;-1 -1;0 -1]; % nodes  Lshaped
 elem = [1 2 3;4 2 3;3 4 5;6 5 4;2 7 4;8 4 7];    % elements
 node = [0 1;0 0.5;0.5 1;0.5 0.5;1 1;1 0.5;1 0;0.5 0;0 0]; % nodes   on the square domain
elem = [1 2 3;4 2 3;3 4 5;6 5 4;4 8 6;7 8 6;2 9 4;8 4 9];    % elements
elem = fixorientation(node,elem);   % counter-clockwise oritentation
elem = label(node,elem);   % label the mesh by the longest edge rule
% showmesh(node,elem);                            % plot mesh
% findelem(nodeH,elemH);                            % plot element indices
% findnode(nodeH);                             % plot node indices
tic
%%  Get a fine mesh by uniform bisection
for k = 1:5
    [node,elem] = uniformrefine(node,elem,[],[]);
end

% for k=1:2
%     [nodeH,elemH,HB,bdFace] = uniformbisect3(nodeH,elemH,HB,bdFace);%showmesh3(node,elem,[130,28],0.35);
% end
%         showmesh3(nodeh,elemh)
% pde = Maxwelldata2;

% bdFlag = setboundary3(node,elem,'Neumann');
bdFlag = setboundary(elem,1);
% bdFlag = setboundary3(node,elem,'Dirichlet','all');   % boundary faces
[elem,bdFlag] = sortelem2(elem,bdFlag);
%not suited for the longest-edge rule
%% Construct Data Structure
% [Dlambda,area] = gradbasis(node,elem);
% [elem2edge,edge,elem2edgeSign] = dofedge(elem);
T= auxstructure(elem(:,1:3)); %elem2edge=T.elem2edge;
edge=T.edge;edge2elem=T.edge2elem; elem2edge = T.elem2edge;
% [dudxc1,dudyc2,jumpdudxc1,jumpdudxc2,jumpdudyc1,jumpdudyc2,area,Dlambda] = gradientconstant(node,elem,u,neighbor);
ve1 = node(elem(:,3),:)-node(elem(:,2),:);
ve2 = node(elem(:,1),:)-node(elem(:,3),:);
ve3 = node(elem(:,2),:)-node(elem(:,1),:);
% elem2edgel=[sum(ve1.*ve1,2) sum(ve2.*ve2,2) sum(ve3.*ve3,2)];
edgel=sqrt(sum((node(edge(:,1),:)-node(edge(:,2),:)).^2,2));
area = 0.5*(-ve3(:,1).*ve2(:,2) + ve3(:,2).*ve2(:,1));
    Dlambda1 = [-ve1(:,2)./(2*area), ve1(:,1)./(2*area)];% [Dlambda1/Dx, Dlambda1/Dy]
    Dlambda2 = [-ve2(:,2)./(2*area), ve2(:,1)./(2*area)];
    Dlambda3 = [-ve3(:,2)./(2*area), ve3(:,1)./(2*area)];
elemsign=sign(area); area=abs(area);
normal(:,:,1)  = [ve1(:,2).*elemsign, -ve1(:,1).*elemsign];% [Dlambda1/Dx, Dlambda1/Dy]
normal(:,:,2)  = [ve2(:,2).*elemsign, -ve2(:,1).*elemsign];
normal(:,:,3)  = [ve3(:,2).*elemsign, -ve3(:,1).*elemsign];
% normal(:,:,1) = -Dlambda1./repmat(sqrt(Dlambda1(:,1).^2 + ...
%                                              Dlambda1(:,2).^2),1,2);%第一条边的单位法向=[ normal(:,1,1), normal(:,2,1) ] 适用于非逆时针排序情况
% normal(:,:,2) = -Dlambda2./repmat(sqrt(Dlambda2(:,1).^2 + ...
%                                              Dlambda2(:,2).^2),1,2);
% normal(:,:,3) = -Dlambda3./repmat(sqrt(Dlambda3(:,1).^2 + ...
%                                              Dlambda3(:,2).^2),1,2);
    % dir=accumarray([face2elem(:,1),face2elem(:,3)],1,[NT,4]);%取边的第一个单元的外法向基函数（囊括了边界边）
% inedgejudge=face2elem(:,1)~=face2elem(:,2);%内部边逻辑向量
% dir=dir+accumarray([face2elem(inedgejudge,2),face2elem(inedgejudge,4)],-1,[NT,4]);
% td=[ones(NT,18),repmat(dir,1,3)]; d=zeros(NT,30,30);
% for i=1:30;
%     for j=1:30
%         d(:,i,j)=td(:,i).*td(:,j);
%     end
% end
B=sym('[B11 B12;B21 B22]');
%det(B)=-volume*(-1).^idx, %volume*(-1).^idx是用12 13 14 vector计算的行列式,
%det(B):41 42 43  vector计算的行列式, 两者正好相差-1
Dl=sym('[D1l1 D2l1;D1l2 D2l2]');
u=Dl.'*basis;%u(1:3,1:30)
c1u=curlbasis;%c1u(1:3,1:30), not divide by |B|=|A31 A32|
c2u=sym(zeros(2,12)); 
for i=1:12;
    Dcurlbasis=[diff(curlbasis(:,i),'lam1') diff(curlbasis(:,i),'lam2')];
    Dxcurlu=Dcurlbasis*Dl;%not divide by |B|
    c2u(:,i)=[Dxcurlu(2);-Dxcurlu(1)];
end
B1=node(elem(:,1),:)-node(elem(:,3),:);  B2=node(elem(:,2),:)-node(elem(:,3),:); 
% B11=B1(:,1); B21=B1(:,2); B12=B2(:,1); B22=B2(:,2); 
Dl1=Dlambda1; Dl2=Dlambda2; 
% D1l1=Dl1(:,1); D2l1=Dl1(:,2); D1l2=Dl2(:,1); D2l2=Dl2(:,2);
%  T= auxstructure(elem(:,1:3)); elem2edge=T.elem2edge;
%  edge=T.edge;edge2elem=T.edge2elem; 
NT = size(elem,1);  %N = size(node,1);
NE = size(edge,1); Ndof = 3*(NE + NT);
elem2dof = [elem2edge elem2edge+NE elem2edge+2*NE];
elem2dof =[double(elem2dof) (1:NT)'+3*NE (1:NT)'+NT+3*NE (1:NT)'+2*NT+3*NE];
    
[lambda,weight] = quadpts(4);lambda1=lambda(:,1)';lambda2=lambda(:,2)';lambda3=lambda(:,3)'; 
lam1=repmat(lambda1,NT,1); lam2=repmat(lambda2,NT,1); 
weight=weight(:)';
nQuadT=size(weight,2); %x0=(node(elem(:,1),1)+node(elem(:,2),1)+node(elem(:,3),1)+node(elem(:,4),1))/4; 
% y0=(node(elem(:,1),2)+node(elem(:,2),2)+node(elem(:,3),2)+node(elem(:,4),2))/4;
% z0=(node(elem(:,1),3)+node(elem(:,2),3)+node(elem(:,3),3)+node(elem(:,4),3))/4;
% inv([16 x y;x 16 z;y z 14]-eye(3))
x=kron(node(elem(:,1),1),lambda1)+kron(node(elem(:,2),1),lambda2)+kron(node(elem(:,3),1),lambda3);
y=kron(node(elem(:,1),2),lambda1)+kron(node(elem(:,2),2),lambda2)+kron(node(elem(:,3),2),lambda3);
% nv11=(z^2 - 195)/(13*x^2 - 2*x*y*z + 15*y^2 + 15*z^2 - 2925);
% nv12=(13*x - y*z)/(13*x^2 - 2*x*y*z + 15*y^2 + 15*z^2 - 2925);
% nv13=(15*y - x*z)/(13*x^2 - 2*x*y*z + 15*y^2 + 15*z^2 - 2925);
% nv21=nv12; nv22=(y^2 - 195)/(13*x^2 - 2*x*y*z + 15*y^2 + 15*z^2 - 2925);
% nv23=(15*z - x*y)/(13*x^2 - 2*x*y*z + 15*y^2 + 15*z^2 - 2925);
% nv31=nv13; nv32=nv23; nv33=(x^2 - 225)/(13*x^2 - 2*x*y*z + 15*y^2 + 15*z^2 - 2925);
nd=16;%n=8+xt-yt;
% y0=repmat(y0,1,nQuadT); x0=repmat(x0,1,nQuadT); z0=repmat(z0,1,nQuadT);
% xt=xt-x0; yt=yt-y0; zt=zt-z0;
B11=repmat(B1(:,1),1,nQuadT); B21=repmat(B1(:,2),1,nQuadT); B12=repmat(B2(:,1),1,nQuadT); 
B22=repmat(B2(:,2),1,nQuadT); 
D1l1=repmat(Dl1(:,1),1,nQuadT); D2l1=repmat(Dl1(:,2),1,nQuadT); 
D1l2=repmat(Dl2(:,1),1,nQuadT); D2l2=repmat(Dl2(:,2),1,nQuadT);  
% lambda=zeros(NT,nQuad,3);
weight=repmat(weight,NT,1); %size(wd)
phi=zeros(NT,nQuadT,12,2); curlphi=phi; curl2phi=phi;  %bt=zeros(NT,nQuadT,3);
for j=1:12
    phi(:,:,j,1)=eval(u(1,j)); phi(:,:,j,2)=eval(u(2,j));
    curlphi(:,:,j)=eval(c1u(j)); 
    curl2phi(:,:,j,1)=eval(c2u(1,j)); curl2phi(:,:,j,2)=eval(c2u(2,j));
end
%% Triangulation output
% T = struct('edge',edge,'face',face,'face2edge',face2edge, ...
%            'elem',elem,'bdFlag',bdFlag);
ii = zeros(144*NT,1); jj = ii; iit = zeros(78*NT,1); jjt = iit; 
index = 0; indext = 0; sA=jj; sAt=sA; sB=iit; sD=sB; sDt=sB; sAc=sB;
% [lambda,w] = quadpts3(3); % quadrature order is 3
% nQuad = size(lambda,1);
for i = 1:12
    for j = 1:12
        ii(index+1:index+NT) = double(elem2dof(:,i));
        jj(index+1:index+NT) = double(elem2dof(:,j));
        %         i1 = locBasesIdx(i,1); i2 = locBasesIdx(i,2); i3 = locBasesIdx(i,3);
        t=(curl2phi(:,:,i,1).*phi(:,:,j,1)+curl2phi(:,:,i,2).*phi(:,:,j,2));
        Aij=sum(t.*weight,2)/2.*elemsign; Atij=sum((nd-1).^(-1).*t.*weight,2)/2.*elemsign;%.*d(:,i,j)
        sA(index+1:index+NT) = Aij; sAt(index+1:index+NT) = Atij; index = index + NT;
        if j>=i
            iit(indext+1:indext+NT) = double(elem2dof(:,i));
            jjt(indext+1:indext+NT) = double(elem2dof(:,j));
            Bij=sum((curlphi(:,:,i,1).*curlphi(:,:,j,1)+curlphi(:,:,i,2).*curlphi(:,:,j,2)).*weight,2)/4./area;
            %          A=A+sparse(elem2dof(:,i),elem2dof(:,j), sum(curlphi(:,:,i).*curlphi(:,:,j).*weight,2).*volume,Ndof,Ndof);
            t=(phi(:,:,i,1).*phi(:,:,j,1)+phi(:,:,i,2).*phi(:,:,j,2));
            Dij=sum(t.*weight,2).*area; Dtij=sum((nd-1).^(-1).*t.*weight,2).*area;
            Acij=sum((nd-1).^(-1).*(curl2phi(:,:,i,1).*curl2phi(:,:,j,1)+curl2phi(:,:,i,2).*curl2phi(:,:,j,2)).*weight,2)/4./area;
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

nQuad=3; QuadRule = gauleg(0,1,nQuad); lambda1=QuadRule.x; w=QuadRule.w;%或者  [x,w]= legs(nQuad);  x=0.5*x+0.5; w=0.5*w;  
lambda1=lambda1(:)';lambda2=1-lambda1;
lambda1=repmat(lambda1,NE,1); lambda2=repmat(lambda2,NE,1); 
lam1=zeros(size(lambda1)); lam2=lam1; 
w=w(:)'; nQuad=size(w,2); w=repmat(w,NE,1); %x0=x0';%x=repmat(x',NT,1);
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
clear phi curlphi;phi{1}=zeros(NE,nQuad,12,2); phi{2}=phi{1};  curlphi{1}=zeros(NE,nQuad,12); curlphi{2}=curlphi{1};curl2phi=phi;
for k=1:2
    %     xtt=(node(elem(face2elem(:,k),1),1)+node(elem(face2elem(:,k),2),1)+node(elem(face2elem(:,k),3),1)+node(elem(face2elem(:,k),4),1))/4;
    %     ytt=(node(elem(face2elem(:,k),1),2)+node(elem(face2elem(:,k),2),2)+node(elem(face2elem(:,k),3),2)+node(elem(face2elem(:,k),4),2))/4;
    %     ztt=(node(elem(face2elem(:,k),1),3)+node(elem(face2elem(:,k),2),3)+node(elem(face2elem(:,k),3),3)+node(elem(face2elem(:,k),4),3))/4;
    t=edge2elem(:,k);
    B1=node(elem(t,1),:)-node(elem(t,3),:);  B2=node(elem(t,2),:)-node(elem(t,3),:); 
    B11=repmat(B1(:,1),1,nQuad); B21=repmat(B1(:,2),1,nQuad); B12=repmat(B2(:,1),1,nQuad);
    B22=repmat(B2(:,2),1,nQuad); 
    Dl1=Dlambda1(t,:); Dl2=Dlambda2(t,:); 
    D1l1=repmat(Dl1(:,1),1,nQuad); D2l1=repmat(Dl1(:,2),1,nQuad); 
    D1l2=repmat(Dl2(:,1),1,nQuad); D2l2=repmat(Dl2(:,2),1,nQuad); 
    for i=1:3
        t=edge2elem(:,k+2)==i;
        switch i
            case 1
                lam1(t,:)=0; lam2(t,:)=lambda1(t,:); 
            case 2
                lam2(t,:)=0; lam1(t,:)=lambda1(t,:); 
            case 3
                lam1(t,:)=lambda1(t,:); lam2(t,:)=lambda2(t,:);
        end
    end
    % xtt=xx-repmat(xtt,1,nQuad); ytt=yy-repmat(ytt,1,nQuad); ztt=zz-repmat(ztt,1,nQuad);
    %     h=power(volume(face2elem(:,k)),1/3); h=repmat(h,1,nQuad);    x=xtt;y=ytt;z=ztt;
    vol=repmat(area(edge2elem(:,k)).*elemsign(edge2elem(:,k)),1,nQuad);
    for j=1:12
        phi{k}(:,:,j,1)=eval(u(1,j)); phi{k}(:,:,j,2)=eval(u(2,j));
        curlphi{k}(:,:,j)=eval(c1u(j))./vol/2;
        curl2phi{k}(:,:,j,1)=eval(c2u(1,j))./vol/2; curl2phi{k}(:,:,j,2)=eval(c2u(2,j))./vol/2;
    end
end
% direction=accumarray([face2elem(:,1),face2elem(:,3)],1,[NT,3]);%取边的第一个单元的外法向基函数（囊括了边界边）
% infacejudge=face2elem(:,1)~=face2elem(:,2);%内部边逻辑向量
% direction=direction+accumarray([face2elem(infacejudge,2),face2elem(infacejudge,4)],-1,[NT,3]);%除边界边外取边的第二个单元的内法向基函数
nx=zeros(NE,nQuad);ny=nx; 
for i=1:3
    t=edge2elem(edge2elem(:,3)==i,1);
    nx(edge2elem(:,3)==i,:)=repmat(normal(t,1,i),1,nQuad); ny(edge2elem(:,3)==i,:)=repmat(normal(t,2,i),1,nQuad);
end
t=edge2elem(:,1)==edge2elem(:,2); %边界边
AE =  sparse(Ndof,Ndof); BE =  AE; AAE =  AE; AE1 =  AE; BE1 =  AE; AAE1 =  AE;
% sAt=zeros(9*dofl^2*Nt,1); mmt=sAt;  nnt=mmt;
for m=1:12%若组装时分成两个单元，组装时间还可以优化
    for n=1:12
%         A1=sum(0.5*curl3phi(:,:,m).*(phi(:,:,n,1).*ny-phi(:,:,n,2).*nx).*w,2).*facel;%单元自身基函数
%         A2=sum(0.5*curl3phi(:,:,n).*(phi(:,:,m,1).*ny-phi(:,:,m,2).*nx).*w,2).*facel;
        A1={0,0;0,0}; BB=A1;A2=A1;% B=A1;C=A1;%curl2uv={0,0};curl3uv={0,0};
%         Amn=sum(0.5*curl3phi(:,:,m).*(phi(:,:,n,1).*ny(:,:,i)-phi(:,:,n,2).*nx(:,:,i)).*w,2).*elem2facel(:,i);%相邻单元之间基函数（方法的优势）
%          Amn=Amn+sum(0.5*curl3phi(:,:,n).*(phi(:,:,m,1).*ny(:,:,i)-phi(:,:,m,2).*nx(:,:,i)).*w,2).*elem2facel(:,i);
        for k1=1:2%1-2号单元
            for k2=1:2
        A1{k1,k2}=sum((-0.5*phi{k1}(:,:,m,1).*ny+0.5*phi{k1}(:,:,m,2).*nx).*curlphi{k2}(:,:,n).*w,2);%.*td(face2elem(:,k1),m).*td(face2elem(:,k2),n);
        A2{k1,k2}=sum((-0.5*curl2phi{k1}(:,:,m,1).*ny+0.5*curl2phi{k1}(:,:,m,2).*nx).*curlphi{k2}(:,:,n).*w,2);
        BB{k1,k2}=sum(curlphi{k1}(:,:,m).*curlphi{k2}(:,:,n).*w,2).*eta1;%单元自身基函数
            end
        end
        AE=AE+sparse(elem2dof(edge2elem(~t,[1 2 1 2]),m),elem2dof(edge2elem(~t,[1 2 2 1]),n),...
            [A1{1,1}(~t); -A1{2,2}(~t);-A1{1,2}(~t); A1{2,1}(~t)],Ndof,Ndof); 
%         AE=AE+sparse(elem2dof(face2elem(t,1),m),elem2dof(face2elem(t,1),n),2*A1{1,1}(t),Ndof,Ndof); 
        AAE=AAE+sparse(elem2dof(edge2elem(~t,[1 2 1 2]),m),elem2dof(edge2elem(~t,[1 2 2 1]),n),...
            [A2{1,1}(~t); -A2{2,2}(~t);-A2{1,2}(~t); A2{2,1}(~t)],Ndof,Ndof); 
%         AAE=AAE+sparse(elem2dof(face2elem(t,1),m),elem2dof(face2elem(t,1),n),2*A2{1,1}(t),Ndof,Ndof); 
        BE=BE+sparse(elem2dof(edge2elem(~t,[1 2 1 2]),m),elem2dof(edge2elem(~t,[1 2 2 1]),n),...
            [BB{1,1}(~t); BB{2,2}(~t);-BB{1,2}(~t); -BB{2,1}(~t)],Ndof,Ndof); 
        AE1=AE1+sparse(elem2dof(edge2elem(t,1),m),elem2dof(edge2elem(t,1),n),2*A1{1,1}(t),Ndof,Ndof); 
        AAE1=AAE1+sparse(elem2dof(edge2elem(t,1),m),elem2dof(edge2elem(t,1),n),2*A2{1,1}(t),Ndof,Ndof); 
        BE1=BE1+sparse(elem2dof(edge2elem(t,1),m),elem2dof(edge2elem(t,1),n),BB{1,1}(t),Ndof,Ndof); 
    end
end
AE=AE1+AE; BE=BE1+BE; AAE=AAE1+AAE;
% AE =  sparse([mm;mmt],[nn;nnt],[sAnt;sAt],Ndof,Ndof); BE =  sparse([mm;mmt],[nn;nnt],[sBnt;sBt],Ndof,Ndof);
% AAE =  sparse([mm;mmt],[nn;nnt],[sAAnt;sAAt],Ndof,Ndof);
%% Boundary conditions
%% Part 1: Find Dirichlet dof and modify the matrix
% Find Dirichlet boundary dof: fixedDof
    % Dirichlet boundary condition only
%     bdFlag = setboundary3(elem,1);
isBdDof = false(NE,1);  
isBdDof(elem2edge(bdFlag==1)) = true;  
t=[~isBdDof;~isBdDof;~isBdDof;true(3*NT,1)]; isBdDof=~t;   
freedof=[~isBdDof;true(Ndof,1);~isBdDof];
    tic;
    AA=[BE/100+Ac+(AAE'+AAE)/15,sparse(Ndof,Ndof);sparse(Ndof,Ndof), speye(Ndof)];
    BB=[At+At'+B+(AE'+AE)/15, -Dt-D; speye(Ndof),sparse(Ndof,Ndof)];
    [eigf,eigv]=eigs(AA([~isBdDof;~isBdDof],[~isBdDof;~isBdDof]),BB([~isBdDof;~isBdDof],[~isBdDof;~isBdDof]),4,1.9^2);
    diag(eigv).^.5
    for i=1:length(eigv)
        indicator(i,:)=abs(  eigf(1:end/2,i)'*(-At(~isBdDof,~isBdDof)-At(~isBdDof,~isBdDof)'+...
            2*eigv(i,i)*(Dt(~isBdDof,~isBdDof)+D(~isBdDof,~isBdDof))-AE(~isBdDof,~isBdDof)/15-AE(~isBdDof,~isBdDof)'/15)*eigf(1:end/2,i)...
        /(eigf(1:end/2,i)'*B(~isBdDof,~isBdDof)*eigf(1:end/2,i))-1  );
    end
    indicator
    toc
    
  