%-------------------------------------------------------------------------
%----- Fichier : deconv_MC
%----- Objet   : Déconvolution par moindres carrés régularisés ou non
%----- Références principales : Idier
%-------------------------------------------------------------------------

close all;
clear all;

L_x = 100;	    % durée d'observation du signal d'entrée
fe = 1;		    % fréquence d'échantillonage
Te = 1/fe;		% période d'échantillonage
N_x = L_x/Te;	% nombre de points du signal x
k_x = 0:N_x;	% index temporel
t_x= k_x*Te;    % base de temps
%**************************************************************************
% 1. Signal x
k1 = 0:N_x/2-1;
t1 = k1*Te;
x1 = t1/(L_x/2);
k2 = N_x/2:N_x;
t2 = k2*Te;
x2 = (-t2 + L_x)/(L_x/2);
x = [x1'; x2'];
subplot(2,3,1)
plot(t_x,x);
xlabel('temps (s)'); ylabel('amplitude'); title('Signal entrée x(t)');


%**************************************************************************
% 2. Réponse impulsionnelle h de filtres linéaires gaussiens

mu_h=15*Te;
% test de plusieurs filtres gaussien ; la largeur de la réponse impulsionnelle va rendre le problème plus ou
% moins difficile
%**************************************************************************
N_sigma=30;
for sigma=1:N_sigma,
   sigma_h=sigma*Te;
   L_h=30*Te;
   t_h=(0:Te:L_h);
   N_h=length(t_h);
   h=(1/(sigma_h*sqrt(2*pi)))*exp(-(((t_h-mu_h)/(sqrt(2)*sigma_h)).*((t_h-mu_h)/(sqrt(2)*sigma_h))));
   % subplot(2,3,2)
   % plot(t_h,h)
   % xlabel('temps (s)'); ylabel('amplitude'); title('Réponse impulsionnelle h(t)');
   
   %-----------------------------------------------------------------------
   % 3. Construction de la matrice H de convolution (Toeplitz)
   N_x=length(x);
   N_y = N_x + N_h -1;
   
   hcol_1 = zeros(1,N_y);
   hcol_1(1:N_h) = hcol_1(1:N_h) + h;
   hlig_1 = zeros(1,N_x);  
   hlig_1(1,1) = hcol_1(1,1);
   H = toeplitz(hcol_1,hlig_1);
   cond_H(sigma)=cond(H);
end

%**************************************************************************
% Construction d'un seul filtre gaussien
%**************************************************************************
sigma=4;
sigma_h=sigma*Te;
L_h=30*Te;
t_h=(0:Te:L_h);
N_h=length(t_h);
h=(1/(sigma_h*sqrt(2*pi)))*exp(-(((t_h-mu_h)/(sqrt(2)*sigma_h)).*((t_h-mu_h)/(sqrt(2)*sigma_h))));
subplot(2,3,2)
plot(t_h,h)
xlabel('temps (s)'); ylabel('amplitude'); title('Réponse impulsionnelle h(t)');

%--------------------------------------------------------------------------
%-
% 3. Construction de la matrice H de convolution (Toeplitz)

N_y = N_x + N_h -1;


hcol_1 = zeros(1,N_y);
hcol_1(1:N_h) = hcol_1(1:N_h) + h;
hlig_1 = zeros(1,N_x);  
hlig_1(1,1) = hcol_1(1,1);
H = toeplitz(hcol_1,hlig_1);
cond_H(sigma)=cond(H);

%**************************************************************************
%4. Convolution sous forme matricielle y_nb = H x

y_nb=H*x;
k=0:1:N_y-1;
t=k*Te;
subplot(2,3,3)
plot(t,y_nb)
xlabel('temps (s)'); ylabel('amplitude'); title('Signal sortie non bruit? y_{nb}(t)');

%**************************************************************************
% 5. Simulation de l'opération de quantification ou bruit additif gaussien

y = round(y_nb*2^10)/2^10;

% RSB=20;
% Py_nb=y_nb'*y_nb/N_y;
% sigma_w=sqrt(Py_nb/(10^(RSB/10)))
% w=sigma_w*randn(N_y,1);
% y=y_nb + w;

% représentation temporelle du signal de sortie bruit?
subplot(2,3,6)
stairs(t,y);
xlabel('temps (s)'); ylabel('amplitude'); title('Signal sortie bruit? y(t)');

%**************************************************************************
% 6. Représentation du bruit (de quantification ou gaussien)
w = y - y_nb;
RSB_y=10*log10((y'*y/N_y)/var(w))
subplot(2,3,5)
plot(t,w); 
xlabel('temps (s)'); ylabel('amplitude'); title('Bruit w(t)');

%**************************************************************************
% 7. Déconvolution l2 par moindres carrés

%--------------------------------------------------------------------------
x_rec= inv(H'*H)*H'*y;
%--------------------------------------------------------------------------

subplot(2,3,4)
plot(t_x,x_rec)
xlabel('temps (s)'); ylabel('amplitude'); title('Signal d''entrée reconstruit par MC');
	 

%**************************************************************************
% 8. Déconvolution par moindres carrés régularisés l2, pénalisation des
% différences premières

% Construction de la matrice D1 de différentiation

dcol_1 = zeros(1,N_x);
dcol_1(1:2) = dcol_1(1:2) + [1 -1];
dlig_1 = zeros(1,N_x);  
dlig_1(1,1) = dcol_1(1,1);
D1 = toeplitz(dcol_1,dlig_1);

% Intervalle de variation du coefficient de regularisation
min_alpha=-7;
pas_alpha=0.1;
max_alpha=+2;
i_alpha=0;

% choix "optimal" pour alpha = 0.0013, ie var_alpha=-2.9
%var_alpha=-2.9;

% ligne suivante (ainsi que "end") ? décommenter pour faire varier alpha

for var_alpha=min_alpha:pas_alpha:max_alpha

alpha=10^var_alpha;
i_alpha=i_alpha+1;
%--------------------------------------------------------------------------
x_rec_l2 (:,i_alpha)= inv(H'*H + alpha * D1'*D1)*H'*y;
%--------------------------------------------------------------------------
err_rec(:,i_alpha)=x-x_rec_l2(:,i_alpha);

Werr_rec(i_alpha)=err_rec(:,i_alpha)'*err_rec(:,i_alpha);
F(i_alpha)=norm(D1*x_rec_l2 (:,i_alpha),2)^2;
G(i_alpha)=norm(y-H*x_rec_l2 (:,i_alpha),2)^2;
end

% trac? du signal reconstruit superpos? au signal recherch?
figure
plot(t_x,x)
hold on
plot(t_x,x_rec_l2,'r')
xlabel('temps (s)'); ylabel('amplitude'); title('Signal d''entrée reconstruit MC régularisés');
hold off

%**************************************************************************
% énergie de l'erreur de reconstruction en fonction du coefficient de régularisation
% ? décommenter lorsque l'on fait varier alpha
figure
var_alpha=min_alpha:pas_alpha:max_alpha;
plot(var_alpha,10*log10(Werr_rec))
xlabel('log10(\alpha)'); ylabel('énergie erreur de reconstruction (dB)'); title('Erreur de reconstruction en fonction de \alpha');
%**************************************************************************
% trac? du signal reconstruit superpos? au signal recherch? pour alpha
% optimal (choix "optimal" pour alpha = 0.0013, ie var_alpha=-2.9)
[W_err_rec_min,i_alpha_opt] = min(Werr_rec);
alpha_opt=10^(min_alpha+pas_alpha*(i_alpha_opt-1))
figure
plot(t_x,x)
hold on
plot(t_x,x_rec_l2(:,i_alpha_opt),'r')
xlabel('temps (s)'); ylabel('amplitude'); title('Signal d''entrée reconstruit régularis? pour alpha optimal');
hold off

%**************************************************************************
% comparaison des conditionnements l2
figure
plot(log10(svd(H'*H)))
hold on
plot(log10(svd(H'*H + alpha * D1'*D1)),'g')
xlabel('index'); ylabel('log10(valeurs singulières)'); 
title('valeurs singulières')
legend(' H^TH','H^TH + \alpha D_1^TD_1');
hold off
display('conditionnement')
cond(H'*H)
max(svd(H'*H))/min(svd(H'*H))
cond(H'*H + alpha * D1'*D1)
%**************************************************************************
% Courbe en L
figure
plot(log10(G),log10(F),'.')
xlabel('log10(critére moindres carrés)'); ylabel('log10(terme de pénalit?)'); title('Courbe en L');



%**************************************************************************
% conditionnement en fonction de la largeur du filtre

figure
sigma=1:N_sigma;
stem(sigma,log10(cond_H))
grid
xlabel('sigma_h'); ylabel('log10(conditionnement de H)'); title('Conditionnement en fonction de la largeur de la réponse impulsionnelle');

