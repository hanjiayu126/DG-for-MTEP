function [node,elem,bdEdge,HB,belong] = uniformrefine(node,elem,bdEdge,HB,belong)
% UNIFORMREFINE divide each traingle into 4 similar triangles. 
%
% [node,elem,bdEdge] = uniformbisect(node,elem,bdEdge)
% [node,elem] = uniformbisect(node,elem)
%
% bdEdge in the input and output is optional and may be omitted.
%--------------------------------------------------------------------------
% Copyright (C) 2008 Long Chen. See COPYRIGHT.txt for details.
%--------------------------------------------------------------------------

%------------------- Construct data structure -----------------------------
totalEdge = sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2);
[edge, i2, j] = unique(totalEdge,'rows','legacy');%
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
elem2edge = reshape(j,NT,3);%elem2edge(i,j):the edge number for the jth edge of the elem i
%------------------- Add new nodes --------------------------
node(N+1:N+NE,:) = (node(edge(:,1),:)+node(edge(:,2),:))/2;
% (i>=N+1) node(i,j) :the added node i is the middle node of the original
% edges
HB(:,1) = (N+1:N+NE)'; HB(:,2) = edge(:,1); HB(:,3) = edge(:,2);
edge2newNode = int32((N+1:N+NE)'); %edge2newNode = int32(edge2newNode);
%--------------------------------------------------------------------------
% refine each triangle into four triangles by regular refinement（相对位置）
%                3
%        /    elem3     \
%edge2  5        -       4      edge1
%      /elem1 \elem4 / elem2\
%     1 -        6 -          2
%-------       edge3----------------------------------------------------------------
t=1:NT;%update elem
p(t,1:3) = elem(t,1:3);%
p(t,4:6) = edge2newNode(elem2edge(t,1:3));%original elem num NT is not updated.length(edge2newNode)=NE
%第1:NE条边中点对应1:NE个new node（编号对应），此处生成节点456所在节点编号（更新后的）
elem(t,:) = [p(t,1), p(t,6), p(t,5)];
%单元1:NT的节点1对应新单元1:NT的第一个点, 新单元1:NT的节点2 3对应节点6 5（对应单元1:NT的边3 2中点）
elem(NT+1:2*NT,:) = [p(t,6), p(t,2), p(t,4)];
elem(2*NT+1:3*NT,:) = [p(t,5), p(t,4), p(t,3)];
elem(3*NT+1:4*NT,:) = [p(t,4), p(t,5), p(t,6)];
%the variable belong is added by han.
if nargin==5
belong(1:4*NT)=belong([t,t,t,t]);
end
%-------------------- Update boundary edges-------------------------------
if ~isempty(bdEdge)%if (nargin==3) && (~isempty(bdEdge))   注意这里被我改了
    bdEdge(NT+1:2*NT,[1 3]) = bdEdge(t,[1 3]); 
    bdEdge(2*NT+1:3*NT,[1 2]) = bdEdge(t,[1 2]); 
    bdEdge(3*NT+1:4*NT,1) = 0;
    bdEdge(t,1) = 0;%bdEdge(1:NT,2:3)不变
else
    bdEdge = [];
end