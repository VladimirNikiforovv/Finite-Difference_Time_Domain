function ep = epsilonCUB(yy)

        if yy > 0 && yy < 0.1
          ep = 81;%2.71 1.97 1200
        else
          ep = 1;
        end


