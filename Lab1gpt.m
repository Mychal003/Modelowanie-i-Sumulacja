%% Laboratorium - Przygotowanie modelu termicznego

clc; clear all; close all;

%% 1. Wyznaczenie parametrów modelu

% Wartości nominalne
TzewN = -20;  % Temperatura zewnętrzna [°C]
Tl_nominal = 20; % Nominalna temperatura w lewym pokoju [°C]
Tp_nominal = 15; % Nominalna temperatura w prawym pokoju [°C]

a = 2; % Współczynnik przenikania ciepła [W/°C]
B = 5; % Grubość ściany działowej [m]
Pgn = 10000; % Moc grzałki [W]

% Wymiary pokoi
x = (50 / B * 2 + 5) / 3;
y = (x - 5) / 2;

Vp = B * x * 3; % Objętość prawego pokoju [m^3]
Vl = B * y * 3; % Objętość lewego pokoju [m^3]

% Parametry powietrza
Cp = 1000; % Ciepło własciwe powietrza [J/(kg*K)]
rop = 1.2; % Gęstość powietrza [kg/m^3]

% Pojemności cieplne
Cvp = Cp * rop * Vp; % Pojemność cieplna prawego pokoju [J/°C]
Cvl = Cp * rop * Vl; % Pojemność cieplna lewego pokoju [J/°C]

%% 2. Obliczenie punktu pracy (punkt równowagi)

% Przewodności cieplne między pomieszczeniami
Ksp = Pgn / (a * (Tl_nominal - TzewN) + (Tp_nominal - TzewN));
Ksl = a * Ksp;
Ksw = Ksp * (Tp_nominal - TzewN) / (Tl_nominal - Tp_nominal); % Przewodność cieplna między pokojami




%--------------------



% Punkt równowagi dla Tl i Tp
Pg0 =Pgn;
Tzew0 = TzewN;

czasskok = 500;

Tl_eq = Tl_nominal;
Tp_eq = Tp_nominal;

%% 3. Sprawdzenie poprawności obliczeń (punkt pracy = wartości nominalne)

dPg = 1;
dTzew =0;
%Tl_eq = (Pg0 + Ksl * Tzew0 + Ksw * Tp_eq) / (Ksl + Ksw)
%Tp_eq = (Ksw * Tl_eq + Ksp * Tzew0) / (Ksw + Ksp)

% Porównanie wartości obliczonych z nominalnymi
disp('Obliczone wartości punktu pracy:');
disp(['Temperatura w lewym pokoju (Tl_eq) = ', num2str(Tl_eq), '°C']);
disp(['Temperatura w prawym pokoju (Tp_eq) = ', num2str(Tp_eq), '°C']);

disp('Nominalne wartości:');
disp(['Temperatura nominalna w lewym pokoju (Tl_nominal) = ', num2str(Tl_nominal), '°C']);
disp(['Temperatura nominalna w prawym pokoju (Tp_nominal) = ', num2str(Tp_nominal), '°C']);

% Sprawdzenie czy obliczone wartości są równe nominalnym
if abs(Tl_eq - Tl_nominal) < 1e-3 && abs(Tp_eq - Tp_nominal) < 1e-3
    disp('Punkt pracy zgadza się z wartościami nominalnymi.');
else
    disp('Punkt pracy NIE zgadza się z wartościami nominalnymi.');
end

%% 4. (Opcjonalnie) Zastosowanie równań stanu

% Macierze równań stanu (przykład z dwoma temperaturami)
A = [-Ksl/Cvl, Ksw/Cvl;
      Ksw/Cvp, -Ksp/Cvp];

B = [1/Cvl; 0];  % Wpływ mocy grzałki na temperaturę Tl
C = eye(2);      % Wyjściem są oba stany (Tl, Tp)
D = [0; 0];      % Brak bezpośredniego wpływu wejść na wyjścia

% Wyświetlenie macierzy równań stanu
disp('Macierz A:');
disp(A);

disp('Macierz B:');
disp(B);

disp('Macierz C:');
disp(C);

disp('Macierz D:');
disp(D);

% Macierze do iteracji (pierwszy etap)
tab_Tzew = [TzewN, TzewN + 10, TzewN - 10];
tab_Tzew1 = [0, -10, 0];% Zmiana temperatury zewnętrznej
tab_Pg = [1, 1, 0.8];             % Zmniejszenie mocy grzałki
tab_color = {'r', 'g', 'b', 'c', 'm', 'k', [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.6350 0.0780 0.1840]};                % Kolory do wykresów

% Wykresy dla pierwszego etapu
figure(1); hold on; grid on; title('Temperatura w lewym pokoju (Tl) - 1 etap');
xlabel('Czas [s]'); ylabel('Temperatura [°C]');

figure(2); hold on; grid on; title('Temperatura w prawym pokoju (Tp) - 1 etap');
xlabel('Czas [s]'); ylabel('Temperatura [°C]');

% Wykresy dla drugiego etapu
figure(3); hold on; grid on; title('Temperatura w lewym pokoju (Tl) - 2 etap');
xlabel('Czas [s]'); ylabel('Temperatura [°C]');

figure(4); hold on; grid on; title('Temperatura w prawym pokoju (Tp) - 2 etap');
xlabel('Czas [s]'); ylabel('Temperatura [°C]');
l=1;
%% Iteracja przez przypadki
for i = 1:3
    Tzew0 = tab_Tzew(i); 
    %Pg0 = Pgn;
    Tl_eq = Tl_nominal;
    Tp_eq = Tp_nominal;
    dPg = 1;
    dTzew = 0;
    out1 = sim("untitled.slx", 'StopTime', '5000'); 
    figure(1);
    plot(out1.tout, out1.Tl_eq, 'Color', tab_color{i}, 'DisplayName', ...
        ['Tzew=', num2str(Tzew0), ', Pg=', num2str(Pg0)]);
    figure(2);
    plot(out1.tout, out1.Tp_eq, 'Color', tab_color{i}, 'DisplayName', ...
        ['Tzew=', num2str(Tzew0), ', Pg=', num2str(Pg0)]);

    for j = 1:3
    dPg = tab_Pg(j);
    dTzew = tab_Tzew1(j); 
    Tl_eq = out1.Tl_eq(end);
    Tp_eq = out1.Tp_eq(end);
    out2 = sim("untitled.slx", 'StopTime', '5000'); 
    figure(3);
    plot(out2.tout, out2.Tl_eq, 'Color', tab_color{l}, 'DisplayName', ...
        ['Tzew=', num2str(Tzew0+dTzew), ', Pg=', num2str(Pg0*dPg)]);

    figure(4);
    plot(out2.tout, out2.Tp_eq, 'Color', tab_color{l}, 'DisplayName', ...
        ['Tzew=', num2str(Tzew0+dTzew), ', Pg=', num2str(Pg0*dPg)]);
    l=l+1;
    end
end

figure(1); legend('show');
figure(2); legend('show');
figure(3); legend('show');
figure(4); legend('show');