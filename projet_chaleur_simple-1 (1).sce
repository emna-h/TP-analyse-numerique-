// =====================================================================
//  TP Analyse numerique - Equation de la chaleur (differences finies)
//  Version "niveau cours" : on n'utilise que ce qui a ete vu
//  (vecteurs, linspace, A\b, A(i,j), size, function, boucles for, plot)
//  Lancer avec : exec("projet_chaleur_simple.sce", -1)
// =====================================================================
clear; clc;

// ---------------------------------------------------------------------
// Donnees numeriques du TP
// ---------------------------------------------------------------------
L = 1;  T0 = 20;  TL = 100;  k = 50;  q = 10000;

// =====================================================================
// PARTIE III : matrice A  (taille (N-1) x (N-1))
// =====================================================================
function A = matrice_chaleur(N)
    n = N - 1;
    A = zeros(n, n);          // matrice nulle de taille n x n
    for i = 1:n
        A(i, i) = -2;         // diagonale principale
        if i < n then
            A(i, i+1) = 1;    // diagonale au-dessus
            A(i+1, i) = 1;    // diagonale en dessous
        end
    end
endfunction

// =====================================================================
// Second membre b (source uniforme q)
// =====================================================================
function b = second_membre(N, L, T0, TL, q, k)
    h = L / N;
    b = -(q * h^2 / k) * ones(N-1, 1);   // terme de la source
    b(1)   = b(1)   - T0;                // condition a gauche
    b(N-1) = b(N-1) - TL;                // condition a droite
endfunction

// Second membre pour une source variable :
// qx = vecteur colonne des valeurs de q(x) aux N-1 noeuds interieurs
function b = second_membre_var(N, L, T0, TL, qx, k)
    h = L / N;
    b = -(h^2 / k) * qx;
    b(1)   = b(1)   - T0;
    b(N-1) = b(N-1) - TL;
endfunction

// =====================================================================
// PARTIE IV : fonction complete de resolution
// =====================================================================
function [x, T] = resolution_chaleur(N, L, T0, TL, q, k)
    h = L / N;                           // 1. pas
    x = (0:N)' * h;                      // 2. maillage (vecteur colonne)
    A = matrice_chaleur(N);              // 3. matrice A
    b = second_membre(N, L, T0, TL, q, k); // 4. second membre
    Tint = A \ b;                        // 5. resolution de A*T = b
    T = [T0; Tint; TL];                  // 6. ajout des conditions limites
endfunction

// =====================================================================
// Solutions exactes
// =====================================================================
// Source uniforme : T = -q/(2k) x^2 + c1 x + T0
function T = T_exact(x, L, T0, TL, q, k)
    c1 = (TL - T0)/L + q*L/(2*k);
    T = -q/(2*k) * x.^2 + c1 * x + T0;
endfunction

// Source q0 (1 + x/L)
function T = T_exact_var(x, L, T0, TL, q0, k)
    c1 = (TL - T0)/L + 2*q0*L/(3*k);
    T = -q0/k * (x.^2/2 + x.^3/(6*L)) + c1 * x + T0;
endfunction

// Source q0 sin(pi x / L)  (experience libre)
function T = T_exact_sin(x, L, T0, TL, q0, k)
    T = q0*L^2/(k*%pi^2) * sin(%pi*x/L) + (TL - T0)*x/L + T0;
endfunction

// =====================================================================
// PARTIE I : questions 1 a 4
// =====================================================================
N = 10;
h = L / N;
disp("Q1 : h = " + string(h));
disp("Q2 : x0 = " + string(0) + ", x1 = " + string(h) + ", x2 = " + string(2*h) + ", x3 = " + string(3*h));
disp("Q3 : 2 temperatures connues (T0 et TN)");
disp("Q4 : " + string(N-1) + " inconnues (N-1)");

// =====================================================================
// Questions 10 et 11 : test de la matrice
// =====================================================================
A5 = matrice_chaleur(5);
disp(A5);
disp(size(A5));                          // doit afficher 4 4
A10 = matrice_chaleur(10);
disp(size(A10));                         // 9 9
A20 = matrice_chaleur(20);
disp(size(A20));                         // 19 19
disp(second_membre(5, L, T0, TL, q, k)); // b pour N = 5 : -28 -8 -8 -108

// =====================================================================
// PARTIE V : premiere simulation (N = 20)
// =====================================================================
[x, T] = resolution_chaleur(20, L, T0, TL, q, k);
xset("window", 1);
plot(x, T, "o-");
xlabel("x (m)");
ylabel("T (degres C)");
title("Temperature numerique, N = 20");

// =====================================================================
// PARTIE VI : comparaison avec la solution analytique
// =====================================================================
xf = linspace(0, L, 200)';               // maillage fin pour la courbe exacte
xset("window", 2);
plot(xf, T_exact(xf, L, T0, TL, q, k), "r-");
plot(x, T, "bo");
xlabel("x (m)");
ylabel("T (degres C)");
legend("Exacte", "Numerique N = 20");

// Question 20 : N = 10, 20, 50, 100
xset("window", 3);
plot(xf, T_exact(xf, L, T0, TL, q, k), "r-");
[x, T] = resolution_chaleur(10, L, T0, TL, q, k);   plot(x, T, "bo");
[x, T] = resolution_chaleur(20, L, T0, TL, q, k);   plot(x, T, "g+");
[x, T] = resolution_chaleur(50, L, T0, TL, q, k);   plot(x, T, "mx");
[x, T] = resolution_chaleur(100, L, T0, TL, q, k);  plot(x, T, "k.");
xlabel("x (m)");
ylabel("T (degres C)");
legend("Exacte", "N=10", "N=20", "N=50", "N=100");

// =====================================================================
// PARTIE VII : erreur maximale (source uniforme)
// =====================================================================
disp("--- Erreur maximale, source uniforme ---");
for N = [10 20 40 80 160]
    [x, T] = resolution_chaleur(N, L, T0, TL, q, k);
    e = abs(T_exact(x, L, T0, TL, q, k) - T);   // erreurs e_i
    EN = max(e);
    disp("N = " + string(N) + "   h = " + string(L/N) + "   E_N = " + string(EN));
end

// =====================================================================
// 8.1 : source variable  q0 (1 + x/L)
// =====================================================================
q0 = 10000;
disp("--- Source variable q0(1+x/L) ---");
Eprec = 0;
for N = [10 20 40 80 160]
    h = L / N;
    x = (0:N)' * h;
    qx = q0 * (1 + x(2:N) / L);                  // source aux noeuds interieurs
    A = matrice_chaleur(N);
    b = second_membre_var(N, L, T0, TL, qx, k);
    T = [T0; A \ b; TL];
    EN = max(abs(T_exact_var(x, L, T0, TL, q0, k) - T));
    if Eprec > 0 then
        p = log(Eprec / EN) / log(2);            // ordre observe
        disp("N = " + string(N) + "   E_N = " + string(EN) + "   p = " + string(p));
    else
        disp("N = " + string(N) + "   E_N = " + string(EN));
    end
    Eprec = EN;
end

// =====================================================================
// PARTIE VIII : etude parametrique
// =====================================================================
// Influence de q (k = 50)
xset("window", 4);
[x, T] = resolution_chaleur(100, L, T0, TL, 0, k);      plot(x, T, "b-");
[x, T] = resolution_chaleur(100, L, T0, TL, 5000, k);   plot(x, T, "g-");
[x, T] = resolution_chaleur(100, L, T0, TL, 10000, k);  plot(x, T, "r-");
[x, T] = resolution_chaleur(100, L, T0, TL, 20000, k);  plot(x, T, "m-");
[x, T] = resolution_chaleur(100, L, T0, TL, 50000, k);  plot(x, T, "k-");
xlabel("x (m)");
ylabel("T (degres C)");
legend("q=0", "q=5000", "q=10000", "q=20000", "q=50000");
title("Influence de q");

// Influence de k (q = 10000)
xset("window", 5);
[x, T] = resolution_chaleur(100, L, T0, TL, q, 10);   plot(x, T, "b-");
[x, T] = resolution_chaleur(100, L, T0, TL, q, 20);   plot(x, T, "g-");
[x, T] = resolution_chaleur(100, L, T0, TL, q, 50);   plot(x, T, "r-");
[x, T] = resolution_chaleur(100, L, T0, TL, q, 100);  plot(x, T, "k-");
xlabel("x (m)");
ylabel("T (degres C)");
legend("k=10", "k=20", "k=50", "k=100");
title("Influence de k");

// Influence du maillage : Tmax, xmax, erreur, temps de calcul
disp("--- Influence du maillage ---");
for N = [10 20 50 100 500]
    tic();
    [x, T] = resolution_chaleur(N, L, T0, TL, q, k);
    temps = toc();
    [Tmax, imax] = max(T);
    erreur = max(abs(T_exact(x, L, T0, TL, q, k) - T));
    disp("N = " + string(N) + "  h = " + string(L/N) + "  Tmax = " + string(Tmax) ..
         + "  xmax = " + string(x(imax)) + "  erreur = " + string(erreur) ..
         + "  temps = " + string(temps) + " s");
end

// =====================================================================
// PARTIE IX : programme de synthese
// =====================================================================
function [x, T, Tmax, xmax] = simulation_chaleur(L, T0, TL, k, q, N)
    // Construction du maillage
    h = L / N;
    x = (0:N)' * h;
    // Construction du systeme
    A = matrice_chaleur(N);
    b = second_membre(N, L, T0, TL, q, k);
    // Resolution
    T = [T0; A \ b; TL];
    // Recherche de Tmax et xmax
    [Tmax, imax] = max(T);
    xmax = x(imax);
    // Trace
    plot(x, T, "b-o");
    xlabel("x (m)");
    ylabel("T (degres C)");
    title("Distribution de temperature");
endfunction

xset("window", 6);
[x, T, Tmax, xmax] = simulation_chaleur(L, T0, TL, k, q, 50);
disp("Tmax = " + string(Tmax) + " C en xmax = " + string(xmax) + " m");

// =====================================================================
// EXPERIENCE LIBRE : source q0 sin(pi x / L)   (ordre 2 visible)
// =====================================================================
disp("--- Experience libre : source q0 sin(pi x/L) ---");
Eprec = 0;
hh = [];  EE = [];                               // pour le trace de convergence
for N = [10 20 40 80 160]
    h = L / N;
    x = (0:N)' * h;
    qx = q0 * sin(%pi * x(2:N) / L);
    A = matrice_chaleur(N);
    b = second_membre_var(N, L, T0, TL, qx, k);
    T = [T0; A \ b; TL];
    EN = max(abs(T_exact_sin(x, L, T0, TL, q0, k) - T));
    hh = [hh, h];
    EE = [EE, EN];
    if Eprec > 0 then
        p = log(Eprec / EN) / log(2);
        disp("N = " + string(N) + "   E_N = " + string(EN) + "   p = " + string(p));
    else
        disp("N = " + string(N) + "   E_N = " + string(EN));
    end
    Eprec = EN;
end

// Trace log-log : on obtient une droite de pente 2
xset("window", 7);
plot(log(hh), log(EE), "o-");
xlabel("log(h)");
ylabel("log(E_N)");
title("Convergence : pente = ordre p");
