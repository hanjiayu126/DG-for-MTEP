function elem = fixorientation(node,elem)
% FIXORIENTATION set all triangles be oritented counter-clockwise
% 
% elem = fixorientation(node,elem)
%
%--------------------------------------------------------------------------
% Copyright (C) 2008 Long Chen. See COPYRIGHT.txt for details.
%--------------------------------------------------------------------------

ve2 = node(elem(:,1),:)-node(elem(:,3),:);
ve3 = node(elem(:,2),:)-node(elem(:,1),:);
area = -ve3(:,1).*ve2(:,2)+ve3(:,2).*ve2(:,1);
ix = (area<0); 
elem(ix,[2 3]) = elem(ix,[3 2]);