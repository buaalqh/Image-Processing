function s=bruitps(e,d)
taille=size(e);
s=e;
for i=1:taille(1)
    for j=1:taille(2)
        p(i,j)=rand; %g¨¦n¨¦rer une variable al¨¦atoire de loi uniforme
        if p(i,j) < d %la probabilit¨¦ p de remplacement des pixels de l¡¯image
        s(i,j)=255;
        end
    end
end