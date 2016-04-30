l = size(s_c,1)

s_mat = zeros(l); 

for i=1:l
    s_mat(i,i) = s_c(i,1); 
end

A = phantom('Modified Shepp-Logan',2000);
h = figure()
subplot(1,2,1), imagesc(A), axis equal tight off
title ('Original matrix')
subplot(1,2,2), imagesc(u_c * s_mat * v_c'), axis equal tight off
title('Sketching approximation')
print -r500
saveas(h,'k_svd.png')