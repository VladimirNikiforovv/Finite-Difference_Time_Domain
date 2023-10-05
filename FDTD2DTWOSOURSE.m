clc
clearvars
cla
close all

%%




x = linspace(-0.1,0.2,300);
y = linspace(-0.1,0.2,300);
T = 600;

dx = x(2)-x(1);
dy = y(2)-y(1);

mu0 = pi * 4e-7;
eps0 = 8.854187817e-12;

c = 1 / sqrt (mu0 * eps0);
% venv = 1 / sqrt(mu0 * eps0 * )
dt = (1 / (sqrt((1 / (dx ^ 2)) + (1 / (dy ^ 2))))) / c;% #* venv;

siz = length(x);

f = 10e9;
phi_0 = 0;
wavelength_m = c / f;
wavelength = wavelength_m / dx;

a = max(x);
b = max(y);

focus = sqrt(a^2 + b^2)/2;

aPix = length(x);
bPix = length(y)/2;
focusPix = sqrt(aPix^2 + bPix^2)/2;

[xelp, yelp] = meshgrid(x, y);

[xx, yy] = meshgrid(x, y);

ugol = linspace(0,2*pi,length(x));

[tetMesh, fiMesh] = meshgrid(ugol, ugol);

for i = 1:siz
    for j = 1:siz
        yelp(i,j) = b * sqrt(1 - (xx(i,j)/a)^2);
    end
end

yEPS= zeros(siz, siz);

for i = 1:siz
    for j = 1:siz
        yEPS(i,j) = b * sqrt(1- (xx(i,j)/a)^2);
    end
end


Hx = zeros(siz, siz - 1, T);
Hy = zeros(siz - 1, siz, T);
Ez = zeros(siz, siz, T);

sigE = zeros(siz, siz);
sigH = zeros(siz, siz);

chxH = zeros(siz, siz);
chxE = zeros(siz, siz);

chyH = zeros(siz, siz);
chyE = zeros(siz, siz);

cEz = zeros(siz, siz);
cHy = zeros(siz, siz);
cHx = zeros(siz, siz);

arrEzProm = zeros(siz, siz, T);

for i = 1:siz
    for j = 1:siz
        chxH(i, j)  = (1 - (dt * sigH(i, j)) / (mu(x(j),y(i)) * mu0 * 2)) / (1 + (dt * sigH(i, j)) / (mu(x(j),y(i)) * mu0 * 2));
        chxE(i, j)  = (dt / (dy * mu(x(j),y(i)) * mu0)) ./ (1 + ((dt * sigH(i, j)) ./ (mu(x(j),y(i)) * mu0 * 2)));

        chyH(i, j) = chxH(i, j);
        chyE(i, j) = (dt / (dx * mu(x(j),y(i)) * mu0)) / (1 + ((dt * sigH(i, j)) / (mu(x(j),y(i)) * mu0 * 2)));


        cEz(i, j) = (1 - (dt * sigE(i, j)) / (epsilon(x(j), y(i)) * eps0 * 2)) / (1 + (dt * sigE(i, j)) / (epsilon(yy(i,j),yEPS(i,j)) * eps0 * 2));
        cHy(i, j) = (dt / (dx * epsilon(yy(i,j),yEPS(i,j)) * eps0)) / (1 + ((dt * sigE(i, j)) / (epsilon(yy(i,j),yEPS(i,j)) * eps0 * 2)));
        cHx(i, j) = (dt / (dy * epsilon(yy(i,j),yEPS(i,j)) * eps0)) / (1 + ((dt * sigE(i, j)) / (epsilon(yy(i,j),yEPS(i,j)) * eps0 * 2)));
    end
end

for t = 1:T-1

    for i = 1:siz
        for j = 1:siz-2
            Hx(i, j + 1, t + 1) = chxH(i, j) * Hx(i, j+1, t) - chxE(i, j) * (Ez(i, j+1, t) - Ez(i, j, t));
        end
    end

    for i = 1:siz-2
        for j = 1:siz
            Hy(i+1, j , t + 1) = chyH(i, j) * Hy(i+1, j , t) + chyE(i, j ) * (Ez(i+1, j , t) - Ez(i, j , t));
        end
    end

    for i = 2:siz-2
        for j = 2:siz-2
            Ez(i, j , t + 1) = cEz(i, j) * Ez(i, j, t) + cHy(i, j ) * (Hy(i+1, j , t + 1 ) - Hy(i, j , t + 1)) - cHx(i, j) * (Hx(i, j+1, t + 1) - Hx(i, j , t + 1));

            if  t < 1100 %i == 3 && (j > 102 && j < siz-302) && t < T/4 %((i >  siz - 8 && i < siz - 3 && j <= 199 && j < 201) && t * dt <= 8*(1/(4*f)) ) %((i> round(siz/2) -20 && i < round(siz/2) + 20 && j <= 4) && t < 16)
                %
                Ez(i, j, t + 1) = Ez(i, j , t+1) + (1/10000)*(initcondGauss(x(i),y(j),0,0.05,dt,0.002,0.002)+...
                    initcondGauss(x(i),y(j),0.1,0.05,dt,0.002,0.002))*...
                    planeWave(0,0,0,0,dt,t+1,1,1,f*27,1000,0);
                %
            end
        end
    end


    Hx(:, 3, t + 1) = Hx(:, 2, t + 1);
    Hx(2, :, t + 1) = Hx(1, :, t + 1);

    Hy(:, 2, t + 1) = Hy(:, 1, t + 1);
    Hy(3, :, t + 1) = Hy(2, :, t + 1);

    Ez(:, 2, t + 1) = Ez(:, 1, t + 1);
    Ez(2, :, t + 1) = Ez(1, :, t + 1);

    Hx(:, siz-1, t + 1) = Hx(:, siz-2, t + 1);
    Hx(siz, :, t + 1) = Hx(siz-1, :, t + 1);

    Hy(:, siz, t + 1) = Hy(:,siz-1, t + 1);
    Hy(siz-1, :, t + 1) = Hy(siz-2, :, t + 1);

    Ez(:, siz-1, t + 1) = Ez(:, siz-2, t + 1);
    Ez(siz-1, :, t + 1) = Ez(siz-2, :, t + 1);
end

% Установка цветовой карты
colormap('gray')

% Создание фигуры
figure(1);

% Создание директории для сохранения файлов
if ~exist('frames', 'dir')
    mkdir('frames');
end

% Цикл по временным шагам
for i = 1:T
    % Отображение текущего состояния
    image(Ez(:,:,i));
    
    % Сохранение текущего изображения в файл
    frame = getframe(gcf);
    imwrite(frame.cdata, sprintf('frames/frame%d.png', i));
    
    % Очистка фигуры
    clf;
end

% Сохранение анимации в файл
filename = 'twoSourse.gif';
for i = 1:T
    % Чтение текущего изображения
    frame = imread(sprintf('frames/frame%d.png', i));
    
    % Преобразование в индексированное изображение
    [imind,cm] = rgb2ind(frame,256);
    
    % Добавление кадра в GIF-файл
    if i == 1
        imwrite(imind,cm,filename,'gif','DelayTime', 0.000001,'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','DelayTime', 0.000001,'WriteMode','append');
    end
end

% Вывод сообщения об успешном сохранении
disp('Анимация успешно записана в файл.');



