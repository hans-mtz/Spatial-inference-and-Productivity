function a = groupcross(g1,g2)

if min(g1) < 1 || min(g2) < 1
    error('Grouping variables must be from 1:G')
end

G2 = max(g2);

a = G2*(g1-1)+g2;
