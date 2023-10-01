function pw = planeWave(x,y,x0,y0,dt,t,kx,ky,omega,A,phi00)

   pw  = A * cos((x-x0)*kx + (y-y0)*ky - omega*t*dt + phi00);
   %pw  = A * exp(