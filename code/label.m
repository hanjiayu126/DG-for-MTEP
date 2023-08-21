function [elem,bdEdge] = label(node,elem,bdEdge)
% LABEL return a tirangulation such that elem(t,1) is opposite to the 
% longest edge of t
%
%    elem = label(node,elem)
%  
%--------------------------------------------------------------------------
% Copyright (C) 2008 Long Chen. See COPYRIGHT.txt for details.
%--------------------------------------------------------------------------

%----------------Compute length of each edge------------------------------
totalEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
ve = node(totalEdge(:,1),:) - node(totalEdge(:,2),:);
edgeLength = reshape(sum(ve.^2,2),size(elem,1),3);
%----------------Switch indices according the edge length-----------------
[temp,I] = max(edgeLength,[],2);
elem((I==2),[1 2 3]) = elem((I==2), [2 3 1]);
elem((I==3),[1 2 3]) = elem((I==3), [3 1 2]);
%------------------- reorder the boundary faces ---------------------------
if nargin==3
	bdEdge((I==2),[1 2 3]) = bdEdge((I==2), [2 3 1]);
	bdEdge((I==3),[1 2 3]) = bdEdge((I==3), [3 1 2]);
else
	bdEdge = [];
end