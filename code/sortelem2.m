function [elem,bdFlag,elem2edge,Dlambda,elem2face] = sortelem2(elem,bdFlag,elem2edge,Dlambda,elem2face)
%% SORTELEM3 sort elem in ascend ordering
%
% [elem,bdFlag] = sortelem3(elem,bdFlag) sorts the elem such that
% elem(t,1)< elem(t,2)< elem(t,3)<elem(t,4). A simple sort(elem,2) cannot
% sort bdFlag.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%%added by Han
if nargin==5
[~,l] = sort(elem,2); locEdge = [1 2; 1 3; 2 3]; NT=size(elem,1);
s=sparse(locEdge(:,1),locEdge(:,2),1:3,3,3); s=s+sparse(locEdge(:,2),locEdge(:,1),1:3,3,3);
elem2edge=elem2edge([(1:NT)' (1:NT)' (1:NT)']+[s(l(:,1)+l(:,2)*3-3)-1 s(l(:,1)+l(:,3)*3-3)-1....
     s(l(:,2)+l(:,3)*3-3)-1]*NT);
% Dlambda1(:,:,1)=reshape(Dlambda(:,1,:),size(elem,1),4);
% Dlambda1(:,:,2)=Dlambda(:,2,:);clear Dlambda
Dlambdax=Dlambda(:,1,:); Dlambday=Dlambda(:,2,:);
Dlambda(:,1,:)=Dlambdax([(1:NT)' (1:NT)' (1:NT)']+[s(l(:,1)+l(:,2)*3-3)-1 s(l(:,1)+l(:,3)*3-3)-1....
     s(l(:,2)+l(:,3)*3-3)-1]*NT);
Dlambda(:,2,:)=Dlambday([(1:NT)' (1:NT)' (1:NT)']+[s(l(:,1)+l(:,2)*3-3)-1 s(l(:,1)+l(:,3)*3-3)-1....
     s(l(:,2)+l(:,3)*3-3)-1]*NT);
% elem2face=[elem2face(l(:,1),:) elem2face(l(:,2),:) elem2face(l(:,3),:) elem2face(l(:,4),:) elem2face(l(:,5),:) elem2face(l(:,6),:)];
end

%% Step 1: elem(:,3) is the biggest one
[tempvar,idx] = max(elem,[],2);  %#ok<*ASGLU>
elem(idx==1,1:3) = elem(idx==1,[2 3 1]);
elem(idx==2,1:3) = elem(idx==2,[3 1 2]);
elem(idx==3,1:3) = elem(idx==3,[2 1 3]);
if exist('bdFlag','var')
    bdFlag(idx==1,1:3) = bdFlag(idx==1,[2  3 1]);
    bdFlag(idx==2,1:3) = bdFlag(idx==2,[3  1 2]);
    bdFlag(idx==3,1:3) = bdFlag(idx==3,[ 2 1 3]);
end
%% Step 2: elem(:,1) is the smallest one
[tempvar,idx] = min(elem(:,1:2),[],2);
% elem(idx==1,1:3) = elem(idx==1,[1 2 3]);
elem(idx==2,1:2) = elem(idx==2,[2 1]);
% elem(idx==3,1:2) = elem(idx==3,[1 2]);
if exist('bdFlag','var')
    bdFlag(idx==2,1:2) = bdFlag(idx==2,[2 1]);
%     bdFlag(idx==3,1:3) = bdFlag(idx==3,[3 1 2]);
end

% %% Step 3: sort elem(:,2)<elem(:3)
% idx = (elem(:,3) < elem(:,2));
% elem(idx,[2 3]) = elem(idx,[3 2]);
% if exist('bdFlag','var')
%     bdFlag(idx,[2 3]) = bdFlag(idx,[3 2]); 
% end

%% Output
if ~exist('bdFlag','var')
    bdFlag = [];
end




 