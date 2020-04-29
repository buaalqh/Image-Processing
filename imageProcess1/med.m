% d¨¦finition de la fonction med (filtrage m¨¦dian)
function s=med(e,voisin)
taille=size(e);
s=e;
Nblim=floor(voisin/2);
for i=1+Nblim:taille(1)-Nblim
    for j=1+Nblim:taille(2)-Nblim
        s(i,j)= median(median(e(i-Nblim:i+Nblim,j-Nblim:j+Nblim)));
    end
end