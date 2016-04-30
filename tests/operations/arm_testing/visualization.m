s_c = importdata('s');
u_c = importdata('u'); 
v_c = importdata('v'); 
l = size(s_c,1)

s_mat = zeros(l); 

for i=1:l
    s_mat(i,i) = s_c(i,1); 
end

A = im2double(rgb2gray(imread('complex2.jpg'))); 
sketch = u_c * s_mat * v_c'; 
h = figure()
subplot(1,2,1), imagesc(A), axis equal tight off
title ('Original matrix')
subplot(1,2,2), imagesc(sketch), axis equal tight off
title('Sketching approximation')