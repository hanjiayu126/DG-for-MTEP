function [elem,bdFace] = label3(node,elem,t,bdFace)
% LABEL3 return a tirangulation such that elem(t,[1 2]) is the longest 
% edge of t.
% 
% Usage
%    elem = label3(node,elem)
%    elem = label3(node,elem,markedElem)
%   [elem,bdFace] = label3(node,elem,t,bdFace)
%
% Input
%    t: selected elements. If it is omitted, t=1:NT.
%   bdFace: when vertices of elements are switched, so is bdFace.


%--------------------------------------------------------------------------
% Copyright (C) 2008 Long Chen. See COPYRIGHT.txt for details.
%--------------------------------------------------------------------------

if (nargin==2)
    t=1:size(elem,1); 
    bdFace = [];
end    
%------------------- Compute edge lengths ---------------------------------
allEdge = [elem(t,[1 2]); elem(t,[1 3]); elem(t,[1 4]); ...
             elem(t,[2 3]); elem(t,[2 4]); elem(t,[3 4])];
allEdgeLength = sum((node(allEdge(:,1),:) - node(allEdge(:,2),:)).^2,2);
elemEdgeLength = reshape(allEdgeLength,length(t),6);
[foo,idx] = max(elemEdgeLength,[],2);
%------------------- reorder the vertices of elem -------------------------
elem(t((idx==2)),1:4) = elem(t((idx==2)),[3 1 2 4]);%点3 1对应第二条边，见allEdge的生成。
elem(t((idx==3)),1:4) = elem(t((idx==3)),[1 4 2 3]);%%%注意生成节点的顺序按照逆时针取
elem(t((idx==4)),1:4) = elem(t((idx==4)),[2 3 1 4]);
elem(t((idx==5)),1:4) = elem(t((idx==5)),[2 4 3 1]);
elem(t((idx==6)),1:4) = elem(t((idx==6)),[4 3 2 1]);
%------------------- reorder the boundary faces ---------------------------
if (nargin==4) && (~isempty(bdFace))%同elem，bdface是个sparse matrix，bdface每行是4个面（或4个点）的编号，实际上存的是点的编号。
    bdFace(t((idx==2)),1:4) = bdFace(t((idx==2)),[3 1 2 4]);
    bdFace(t((idx==3)),1:4) = bdFace(t((idx==3)),[1 4 2 3]);
    bdFace(t((idx==4)),1:4) = bdFace(t((idx==4)),[2 3 1 4]);
    bdFace(t((idx==5)),1:4) = bdFace(t((idx==5)),[2 4 3 1]);
    bdFace(t((idx==6)),1:4) = bdFace(t((idx==6)),[4 3 2 1]);
else
    bdFace = [];
end