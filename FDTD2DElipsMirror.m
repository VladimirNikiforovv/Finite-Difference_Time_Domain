clc
clearvars
cla
close all

%%




x = linspace(-0.1,0.1,300);
y = linspace(-0.1,0.1,300);
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
        yEPS(i,j) = sqrt(0.0001 - (xx(i,j))^2);
    end
end


Hx = zeros(siz, siz, T);
Hy = zeros(siz, siz, T);
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



for i = 1:siz
    for j = 1:siz
        chxH(i, j)  = (1 - (dt * sigH(i, j)) / (mu(x(j),y(i)) * mu0 * 2)) / (1 + (dt * sigH(i, j)) / (mu(x(j),y(i)) * mu0 * 2));
        chxE(i, j)  = (dt / (dy * mu(x(j),y(i)) * mu0)) ./ (1 + ((dt * sigH(i, j)) ./ (mu(x(j),y(i)) * mu0 * 2)));
        
        chyH(i, j) = chxH(i, j);
        chyE(i, j) = (dt / (dx * mu(x(j),y(i)) * mu0)) / (1 + ((dt * sigH(i, j)) / (mu(x(j),y(i)) * mu0 * 2)));
        
        
        cEz(i, j) = (1 - (dt * sigE(i, j)) / (epsilon(yy(i,j),yEPS(i,j)) * eps0 * 2)) / (1 + (dt * sigE(i, j)) / (epsilon(yy(i,j),yEPS(i,j)) * eps0 * 2));
        cHy(i, j) = (dt / (dx * epsilon(yy(i,j),yEPS(i,j)) * eps0)) / (1 + ((dt * sigE(i, j)) / (epsilon(yy(i,j),yEPS(i,j)) * eps0 * 2)));
        cHx(i, j) = (dt / (dy * epsilon(yy(i,j),yEPS(i,j)) * eps0)) / (1 + ((dt * sigE(i, j)) / (epsilon(yy(i,j),yEPS(i,j)) * eps0 * 2)));
    end
end

% проверка корректности задания эпсилона
% arrEPS = zeros(siz, siz);
% 
% for i = 1:siz
%     for j = 1:siz
%         arrEPS(i,j) = epsilon(yy(i,j), yEPS(i,j));
%     end
% end
% 
% mesh(arrEPS);


% Ez(249 + 124, 249, 1) = 2000; 
% for i = 2:siz-1
%     for j = 2:siz-1
%         if (yy(i,j) <= yelp(i,j)) && (yy(i,j) >= -yelp(i,j))
%           Ez(i, j, 1) = initcondGauss(x(i),y(j),0.05,0,dt,0.002,0.002);
%         else
%           Ez(i, j, 1) = 0;
%         end
%     end
% end
% 
% for i = 1:siz-1
%     for j = 1:siz-1
%         if (yy(i,j) <= yelp(i,j)) && (yy(i,j) >= -yelp(i,j)) && i>10 && i < 40 && j > 250 && j < 255 
%           Ez(i, j, 1) = planeWave(x(j),y(j),0,0,dt,1,0,0.5,f*2,1000,phi_0);
%         else
%           Ez(i, j, 1) = 0;
%         end
%     end
% end


% for i = 1:siz-1
%     for j = 1:siz-1
%         if (yy(i,j) <= yelp(i,j)) && (yy(i,j) >= -yelp(i,j)) && (i>180+50 && i < 182+50 && j > 150 - 30 && j <= 150 + 30) 
%           Ez(i, j, 1) = planeWave(0,0,0,0,dt,1,0,0,f*50,2000,phi_0);
%         else
%           Ez(i, j, 1) = 0;
%         end
%     end
% end


for t = 1:T-1 
  for i = 1:siz-1
    for j = 1:siz-1
      if (yy(i,j) <= yelp(i,j)) && (yy(i,j) >= -yelp(i,j))
        Hx(i, j + 1, t + 1) = chxH(i, j) * Hx(i, j+1, t) - chxE(i, j) * (Ez(i, j+1, t) - Ez(i, j, t));

        Hy(i+1, j , t + 1) = chyH(i, j) * Hy(i+1, j , t) + chyE(i, j ) * (Ez(i+1, j , t) - Ez(i, j , t));

        Ez(i, j , t + 1) = cEz(i, j) * Ez(i, j, t) + cHy(i, j ) * (Hy(i+1, j , t + 1 ) - Hy(i, j , t + 1)) - cHx(i, j) * (Hx(i, j+1, t + 1) - Hx(i, j , t + 1));    

        if (i>180+50 && i < 182+50 && j > 150 - 30 && j <= 150 + 30) 
          Ez(i, j, t + 1) = Ez(i, j , t+1) + planeWave(0,0,0,0,dt,t+1,0,0,f*50,1000,phi_0);
        end
        
      else
        Hx(i, j + 1, t + 1) = 0;

        Hy(i+1, j , t + 1) = 0;

        Ez(i, j , t + 1) = 0;
      end
    end
  end
end

% Сохранение анимации
% writerObj = VideoWriter('animation.avi');
% open(writerObj);
% figure(1);
% for i = 1:T
% image(Ez(:,:,i));
% % view(0,90);
% pause(0.00001);
% hold off
% end
% close(writerObj);


figure(1);
for i = 1:T
image(Ez(:,:,i));
% view(0,90);
pause(0.00001);
hold off
end



% figure(2)
% 
% for k = 1:T 
% contour( Ez(:,:,k),4);
% pause(0.0000001);
% hold off
% end
