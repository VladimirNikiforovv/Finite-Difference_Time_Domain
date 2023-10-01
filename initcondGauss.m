function res = initcondGauss(x,y,x0,y0,dt,sigmaX,sigmaY)

   res = (1/(2*pi*sigmaX*sigmaY)) * exp(-(((x-x0)^2/(2*sigmaX^2)) + ((y-y0)^2/(2*sigmaY^2))));