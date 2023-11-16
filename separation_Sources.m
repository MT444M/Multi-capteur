%%   <<<<<<<<<<<<<<<<<<<<<<<<<Multi-Capt>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clc;
close all;
clear all;

% --------------------- SEPARATIN DE SOURCES -----------------------------

%% Partie 1 : Création de signaux et mélange

% Paramètres
fs = 441; % Fréquence d'échantillonnage (en Hz)
t = 0:1/fs:1; % Vecteur de temps de 0 à 1 seconde
% Création du signal sinusoïdal
f1 = 44; % Fréquence du signal sinusoïdal (en Hz)
A1 = 2; % Amplitude
signal_sinusoide = A1 * sin(2 * pi * f1 * t);
signal_sinusoide = signal_sinusoide / mean(abs(signal_sinusoide).^2); % Normalisation

% Création du signal aléatoire
signal_aleatoire = randn(size(t)); % Signal aléatoire gaussien
signal_aleatoire = signal_aleatoire / mean(abs(signal_aleatoire).^2); % Normalisation

% Matrice de mélange
M = [1, 0.5; 0.5, 1]; % Matrice de mélange 
%M inversible car signaux indépendantes

% Mélange des signaux
y = M * [signal_sinusoide; signal_aleatoire];

% Affichage des signaux mélangés
figure;
subplot(4, 1, 1);
plot(t, signal_sinusoide);
title('Signal Sinusoïdal');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(4, 1, 2);
plot(t, signal_aleatoire);
title('Signal Aléatoire');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(4, 1, 3);
plot(t, y(1, :),'r'); % Signal mélangé 1
title('Signaux Mélangés 1');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(4, 1, 4);
plot(t, y(2, :), 'r'); % Signal mélangé 2
title('Signaux Mélangés 2');
xlabel('Temps (s)');
ylabel('Amplitude');


%%  Partie 2 : Séparation des sources
% 1) On connait M
 %z = inv(M)*y ;

%2)       --------Décorrélation des sources sans Bruits--------------
%             ----                ------
% s(t) -------  M ------y(t)-----   B  ----- x(t) Rx 
%             ----                ------
% Calcul de la matrice de covariance Ry
Ry = cov(y');

% Décomposition en valeurs propres de Ry
[Q, D] = eig(Ry);

% Matrice de Décorrélation
B = Q*(D^-1/2)*Q';

x = B*y;

figure()
subplot(2, 1, 1);
plot(t,x(1, :),'r'); % Signal décorrélé 1
title('Signaux décorrélé 1:  Matrix de décorrélation');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(t,x(2, :), 'r'); % Signal décorrélé 2
title('Signaux décorrélé 2: Matrix de décorrélation');
xlabel('Temps (s)');
ylabel('Amplitude');

%%    ------- Séparation des sources: R_x(tau) et Cx ---------------

tau = 0.3; % Décalage temporel en secondes
% Calculez le décalage en nombre d'échantillons
tau_samples = round(tau * fs);

% Créez le signal y(t - tau) en décalant le signal y(t) dans le temps
y_tau = y(:,tau_samples + 1:end); % Supprimez les échantillons décalés

% Calculez les vecteurs moyens des signaux y_t et y_tau
mean_y = mean(y, 2); % La moyenne est prise sur la deuxième dimension
mean_y_tau = mean(y_tau, 2);

% Calculez la covariance entre y_t et y_t_tau
%covariance = sum((y(:,1:end-tau_samples) - mean_y) .* (y_tau - mean_y_tau)) / length(y);

% Calculez la matrice de covariance entre les signaux
covariance_matrix = (1 / (size(y, 2) - 1)) * ((y(:,1:end-tau_samples) - mean_y) * (y_tau - mean_y_tau)');


T_x = covariance_matrix;
T_x = (T_x+T_x')/2;
[Q_Tx, D_Tx] = eig(T_x);

%Separation de sources: Matrice de Covariance
x = B*y;
z = Q_Tx'*x;  % séparation de sources

figure()
subplot(2, 1, 1);
plot(t,z(1, :),'r'); % Signal décorrélé 1
title('Signal décorrélé 1');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(t,z(2, :), 'r'); % Signal décorrélé 2
hold on;
plot(t, -signal_sinusoide)
title('Séparation de source: Matrice de Covariance R_x(30%) ');
legend('signal estimé','signal d''entrée')
xlabel('Temps (s)');
ylabel('Amplitude');

%% ERREUR DE Séparation : EQM
EQM = rmse(z(2, :),-signal_sinusoide)

%% -------------Matrice Globale--------------
%A = M*B*Q_Tx'

%% ------- Séparation des sources:  Cx =  ---------------
x_ = x;
Cum = zeros(2,2,2,2);
for i =1:2
    for j = 1:2
        for k=1:2
            for l = 1:2
                E1 = mean(x_(i,:).*x(j,:).*x_(k,:).*x_(l,:));
                E2 = mean(x_(i,:).*x_(j,:))*mean(x_(k,:).*x_(l,:));
                E3 = mean(x_(i,:).*x_(k,:))*mean(x_(j,:).*x_(l,:));
                E4 = mean(x_(i,:).*x_(l,:))*mean(x_(j,:).*x_(k,:));
                Cum(i,j,l,k)  = E1- E2- E3- E4;
            end
        end
    end
end

Cum1 = Cum(:,:,1,1);
[U,C_s] = eig(Cum1);

T_x_2 = Cum1;
T_x_2 = (T_x_2+T_x_2')/2;
[Q_Tx2, ~] = eig(T_x_2);

z2 = Q_Tx2'*y;  % Q_Tx 
figure()
subplot(2, 1, 1);
plot(z2(2, :),'r'); % Signal décorrélé 1
title('Signal décorrélé 1');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(t,z2(1, :), 'r'); % Signal décorrélé 2
hold on;
plot(t, -signal_sinusoide)
title('Séparation de source: Cumulant ');
legend('signal estimé','signal d''entrée')
xlabel('Temps (s)');
ylabel('Amplitude');

EQM_2_ = rmse(z2(1, :),-signal_sinusoide);
EQM_2 = rmse((1/ max(z2(1, :)))*z2(1, :),-signal_sinusoide)

%% 3)       --------Décorrélation des sources Avec Bruits--------------

%  --------------------- Ajout de bruit sur y ---------------------
% Paramètres
nb_samples = length(t); % Nombre d'échantillons
num_sources = 2;
% Créez les signaux de bruit décorrélés et de même puissance
nb_bruit = 5; % Nombre de signaux de bruit
bruit_power =0.01; % Puissance du bruit (ajustez selon vos besoins)
bruits = sqrt(bruit_power) * randn(nb_bruit, nb_samples); % Génération de bruit gaussien
%bruits = bruits/mean(abs(bruits).^2);
% Mélange des signaux sources avec le bruit
y_2 = zeros(nb_bruit, nb_samples);
for i = 1:num_sources
    y_2(i, :) = M(i, :) * [signal_sinusoide; signal_aleatoire] + bruits(i, :) ;
end
for i = num_sources+1:nb_bruit
    y_2(i, :) = bruits(i, :);
end


% Calcul de la matrice de covariance Ry
Ry_2 = cov(y_2');

% Décomposition en valeurs propres de Ry
[Q, D] = eig(Ry);
[Q_2, D_2] = eig(Ry_2);

d = diag(D_2); % Extraire la diagonale de C
d = flip(d);
figure()
plot(d);
grid on;

% Calcul de la matrice B: B = VDyQs'
Qs = Q_2(:,nb_bruit-1:nb_bruit);
D_2_ = D_2(nb_bruit-1:nb_bruit,nb_bruit-1:nb_bruit); %Dyb(Nsigma+1:Ncap,Nsigma+1:Ncap); %on considere que la partie signal
B_2 = (D_2_^-1/2)*Qs';


%Signaux décor
x_2 = B_2*y_2;

% figure()
% subplot(2, 1, 1);
% plot(t,x_2(1, :)); % Signal décorrélé 1
% title('Signaux décorrélé 1');
% xlabel('Temps (s)');
% ylabel('Amplitude');
% 
% subplot(2, 1, 2);
% plot(t,x_2(2, :), 'r'); % Signal décorrélé 2
% title('Signaux décorrélé 2');
% xlabel('Temps (s)');
% ylabel('Amplitude');

for i =1:2
    for j = 1:2
        for k=1:2
            for l = 1:2
                E1 = mean(x_2(i,:).*x_2(j,:).*x_2(k,:).*x_2(l,:));
                E2 = mean(x_2(i,:).*x_2(j,:))*mean(x_2(k,:).*x_2(l,:));
                E3 = mean(x_2(i,:).*x_2(k,:))*mean(x_2(j,:).*x_2(l,:));
                E4 = mean(x_2(i,:).*x_2(l,:))*mean(x_2(j,:).*x_2(k,:));
                Cum_2(i,j,l,k)  = E1- E2- E3- E4;
            end
        end
    end
end

Cum1_2 = Cum_2(:,:,1,1);
[U_2,C_s_2] = eig(Cum1_2);

T_x_2 = Cum1_2 ;
T_x_2 = (T_x_2+T_x_2')/2;
[Q_Tx2_2, ~] = eig(T_x_2);

z2_2 = Q_Tx2_2'*x_2;  % Q_Tx 
figure()
% subplot(2, 1, 1);
% plot(z2_2(2, :),'r'); % Signal décorrélé 1
% title('Signal décorrélé 1');
% xlabel('Temps (s)');
% ylabel('Amplitude');
% 
% subplot(2, 1, 2);
plot(t,z2_2(1, :)); % Signal décorrélé 2
hold on;
plot(t, -signal_sinusoide)
title('Séparation de source cas BIG DATA: Cumulant ');
legend('signal estimé','signal d''entrée')
xlabel('Temps (s)');
ylabel('Amplitude');

EQM_3 = rmse(z2_2(1, :),-signal_sinusoide)

