function kb=egal(image,N,x)
%% définition de la fonction d’égalisation d’histogramme: ? egal?
taille=size(image);
f=cumsum(N);
%kb是优化后的niveau de gris
%imread进来的数据就是niveau de gris
kb1=255*f/(taille(1)*taille(2));
for i=1:taille(1)
    for j=1:taille(2)
        %matlab 首元素不是0, 1~256
        kb(i,j)=kb1(round(image(i,j)+1));
    end
end