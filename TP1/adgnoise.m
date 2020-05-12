function image = adgnoise(image_nb,RSB)
%add bruit additif gaussien
N=length(image_nb);
Py_nb=image_nb'*image_nb/N;
sigma_w=sqrt(Py_nb/(10^(RSB/10)));
w=sigma_w*randn(N,1);
image=image_nb + w;
end

