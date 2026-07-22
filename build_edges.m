% =========================================================================
%  EDGE TABLE
% =========================================================================
function [edges, elemEdges, edgeSigns] = build_edges(elements)
numElem  = size(elements,1);
edgeMap  = containers.Map('KeyType','char','ValueType','double');
edges    = zeros(3*numElem, 2);
elemEdges = zeros(numElem, 3);
edgeSigns = zeros(numElem, 3);
cnt = 0;
lp  = [1 2; 2 3; 3 1];
for e = 1:numElem
    tri = elements(e,:);
    for k = 1:3
        n1 = tri(lp(k,1));
        n2 = tri(lp(k,2));
        s  = sort([n1 n2]);
        key = sprintf('%d_%d', s(1), s(2));
        if ~isKey(edgeMap, key)
            cnt = cnt + 1;
            edgeMap(key) = cnt;
            edges(cnt,:) = s;
        end
        elemEdges(e,k) = edgeMap(key);
        edgeSigns(e,k) = sign(n2 - n1);   % +1 if local orientation matches global
    end
end
edges = edges(1:cnt, :);
end

