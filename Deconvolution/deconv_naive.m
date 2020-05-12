%-------------------------------------------------------------------------
%----- Fichier : deconv_naive
%----- Objet   : D�convolution naive
%----- R�f�rences principales : Idier, 
%-------------------------------------------------------------------------

close all;
clear all;

L1x = 100;	    % dur�e d'observation du signal d'entr�e (en s)
fe = 1;		    % fr�quence d'�chantillonnage (Hz)
Te = 1/fe;		% p�riode d'�chantillonnage
N1x = L1x/Te;	% nombre de points du signal x
kx = 0:N1x;	    % vecteur des indices temporels  
tx= kx*Te;		% vecteur des instants d'�chantillonnage


% 1. R�ponse impulsionnelle h d'un filtre gaussien
%
mu_h=15*Te;
% la largeur de la r�ponse impulsionnelle va rendre le probl�me plus ou
% moins difficile
sigma_h=5*Te;
tt=(0:Te:30*Te);
h=(1/(sigma_h*sqrt(2*pi)))*exp(-(((tt-mu_h)/(sqrt(2)*sigma_h)).*((tt-mu_h)/(sqrt(2)*sigma_h))));

figure(1);
subplot(2,3,2)
plot(tt,h)
xlabel('temps (s)'); ylabel('amplitude'); title('R�ponse impulsionnelle h(t)');

% 2. Signal x
%
k1 = 0:N1x/2-1;
t1 = k1*Te;
x1 = t1/(L1x/2); 
k2 = N1x/2:N1x;  % modif  N1x-1 => N1x
t2 = k2*Te;
x2 = (-t2 + L1x)/(L1x/2);
x = [x1 x2];
subplot(2,3,1)
plot(tx,x);
xlabel('temps (s)'); ylabel('amplitude'); title('Signal entr�e x(t)');

%3. Convolution
%

y_nb=conv(x,h);
N1=length(y_nb);
k=0:1:N1-1;
t=k*Te;
subplot(2,3,3)
plot(t,y_nb)
xlabel('temps (s)'); ylabel('amplitude'); title('Signal sortie non bruit� y_{nb}(t)');

% 4. Simulation de l'op�ration de quantification
%
y = round(y_nb*2^10)/2^10;  %  q=1/2^10
subplot(2,3,6)
stairs(t,y);
xlabel('temps (s)'); ylabel('amplitude'); title('Signal sortie bruit� y(t)');

% 5. Repr�sentation du bruit de quantification
%
w = y - y_nb;
subplot(2,3,5)
plot(t,w); 
xlabel('temps (s)'); ylabel('amplitude'); title('Bruit w(t)');
RSB=10*log10(y*y'/var(w)) % 



% 6. R�ponse en fr�quence du filtre
%
N=1024;
n = -N/2:N/2-1;
f = n*fe/N;
H = fft(h,N);
figure(2); subplot(2,3,2)
plot(f,fftshift(20*log10(abs(H))));
xlabel('fr�quence (Hz)'); ylabel('dB'); title('dse du filtre H(f)');


% 7. TFD X du signal d'entr�e, TFD Y du signal de sortie, TFD X_rec du signal reconstruit 
%
X = fft(x,N);
subplot(2,3,1)
plot(f,fftshift(20*log10(abs(X))));
xlabel('fr�quence (Hz)'); ylabel('dB'); title('dsp X(f)');

%Y = fft(y,N); % spectre du signal bruite
Y=fft(y_nb,N); % spectre du signal non bruite
subplot(2,3,3)
plot(f,fftshift(20*log10(abs(Y))));
xlabel('fr�quence (Hz)'); ylabel('dB'); title('dsp Y(f)');


% R�ponse en fr�quence du filtre inverse
%
Hinv = 1./H;
hinv = real(ifft(Hinv));
subplot(2,3,5)
plot(f,fftshift(20*log10(abs(Hinv))));
xlabel('fr�quence (Hz)'); ylabel('dB'); title('dse du filtre inverse 1./H');

% R�ponse en fr�quence du bruit
%
W = fft(w,N);
subplot(2,3,6)
plot(f,fftshift(20*log10(abs(W))));
xlabel('fr�quence (Hz)'); ylabel('dB'); title('dsp W(f)');

% calcul du signal reconstruit par filtrage inverse
%
X_rec2 = Y./H;
x_rec2 = real(ifft(X_rec2));

Nrec2=length(x_rec2);
krec2=0:1:Nrec2-1;
trec2=krec2*Te;
figure(1);
subplot(2,3,4)
plot(trec2(1:length(x)),x_rec2(1:length(x)))
xlabel('temps (s)'); ylabel('amplitude'); title('Signal d''entr�e reconstruit par filtrage inverse');



figure(2)
subplot(2,3,3)
plot(f,fftshift(20*log10(abs(Y))));
xlabel('fr�quence (Hz)'); ylabel('dB'); title('dsp Y(f)');


subplot(2,3,1)
plot(f,fftshift(20*log10(abs(X))));
xlabel('fr�quence (Hz)'); ylabel('dB'); title('dsp X(f)');

subplot(2,3,2)
plot(f,fftshift(20*log10(abs(H))));
xlabel('fr�quence (Hz)'); ylabel('dB'); title('dse du filtre H(f)');

subplot(2,3,6)
plot(f,fftshift(20*log10(abs(W))));
xlabel('fr�quence (Hz)'); ylabel('dB'); title('dsp W(f)');

subplot(2,3,5)
plot(f,fftshift(20*log10(abs(Hinv))));
xlabel('fr�quence (Hz)'); ylabel('dB'); title('dse du filtre inverse 1./H');

subplot(2,3,4)
plot(f,fftshift(20*log10(abs(X_rec2))));
xlabel('fr�quence (Hz)'); ylabel('dB'); title('dsp X_{rec2}(f)');

% figure
% subplot(2,3,3)
% pwelch(y,[],[],[],fe,'centered')
% 
% subplot(2,3,1)
% pwelch(x,[],[],[],fe,'centered')
% 
% subplot(2,3,2)
% pwelch(h,[],[],[],fe,'centered')
% 
% subplot(2,3,6)
% pwelch(w,[],[],[],fe,'centered')
% 
% subplot(2,3,5)
% pwelch(hinv,[],[],[],fe,'centered')
% 
% subplot(2,3,4)
% pwelch(x_rec2,[],[],[],fe,'centered')

%R�ponse impulsionnelle du filtre inverse
% hinv = real(ifft(Hinv));
% kh=0:1:N-1;
% th=kh*Te;
% figure; plot(th,hinv);
% xlabel('temps (s)'); ylabel('amplitude'); title('hinv');
% 
% y = real(ifft(Y));
% figure; plot(t,y);
% xlabel('temps (s)'); ylabel('amplitude'); title('Signal y');


% TFD de X � partir de H et Y et repr�sentation temporelle de x
% Xrex = Y./H;
% Xrexvisu = fftshift(abs(Xrex));
% figure; plot(f,Xrexvisu);
% xlabel('fr�quence (Hz)'); ylabel('amplitude'); title('Spectre d''amplitude de Xrex');
% 
% xrex = real(ifft(Xrex));
% figure; plot(t,xrex);
% xlabel('temps (s)'); ylabel('amplitude'); title('Signal xrex');


% Comparaison
% A ce stade, x et xrex sont rigoureusement �gaux : c'est tout � fait normal.
% Voyons maintenant ce qui se passe lorsque le signal y est bruit�.


% 7. Bruit de quantification
% yb = round(y*2^10)/2^10;
% figure; plot(t,yb);
% xlabel('temps (s)'); ylabel('amplitude'); title('Signal y bruit�');
% 
% Yb = fft(yb);
% Ybvisu = fftshift(abs(Yb));
% figure; plot(f,Ybvisu);
% xlabel('fr�quence (Hz)'); ylabel('amplitude'); title('Spectre d''amplitude de y bruit�');

% Xrex = Yb./H;
% Xrexvisu = fftshift(abs(Xrex));
% figure; plot(f,Xrexvisu);
% xlabel('fr�quence (Hz)'); ylabel('amplitude'); title('Spectre d''amplitude de xrex (y bruit�)');
% 
% xrex = real(ifft(Xrex));
% figure; plot(t,xrex); hold on; plot(t,xrex,'r');
% xlabel('temps (s)'); ylabel('amplitude'); title('Signal xrex (y bruit�)');



