function kb=egal(image,N,x)
%% d��finition de la fonction d����galisation d��histogramme: ? egal?
taille=size(image);
f=cumsum(N);
%kb���Ż����niveau de gris
%imread���������ݾ���niveau de gris
kb1=255*f/(taille(1)*taille(2));
for i=1:taille(1)
    for j=1:taille(2)
        %matlab ��Ԫ�ز���0, 1~256
        kb(i,j)=kb1(round(image(i,j)+1));
    end
end