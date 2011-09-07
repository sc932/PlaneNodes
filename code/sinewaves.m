clear
randn('state',sum(100*clock));

relsize = 1000;
Lx = 1.0;
Ly = 1.61803399;
nmax = relsize;
mmax = floor(relsize*sqrt(Ly));

E = zeros(nmax,mmax);
k = 0;
for n = 1:nmax
    for m = 1:mmax
        E(n,m) = (n*n/Lx + m*m/Ly);
        k = k + 1;
        nmtracker(k,:) = [n,m];
    end
end

% Sorts the energies, but keeps track of corresponding n,m
Elist = [];
for i = 1:nmax
    Elist = [Elist,E(i,:)];
end

for i = 1:(nmax*mmax)
    currhigh = 0;
    Ehigh = 0;
    for j = 1:(nmax*mmax - i + 1)
        if (Elist(j) >= Ehigh)
            currhigh = j;
            Ehigh = Elist(j);
        end
    end
    if currhigh > 0
        temp = Elist(nmax*mmax - i + 1);
        Elist(nmax*mmax - i + 1) = Elist(currhigh);
        Elist(currhigh) = temp;L
        temp2 = nmtracker(nmax*mmax - i + 1,:);
        nmtracker(nmax*mmax - i + 1,:) = nmtracker(currhigh,:);
        nmtracker(currhigh,:) = temp2;
    end
end

for i = 1:(nmax*mmax)
    deltaE(i) = sqrt(2*Elist(i))/((Lx+Ly)/2);
end

% plot of E with an error bar of size delta E
errorbar(Elist,deltaE,'o')

% neff finder
killpoint = 0;
for i = 1:(nmax*mmax)
    for j = 1:i
        if((Elist(i) - Elist(j)) < deltaE(i))
            maxawayL(i) = i-j;
            break;
        end
    end
    for j = i:(nmax*mmax)
        if((Elist(j) - Elist(i)) < deltaE(i))
            maxawayU(i) = j-i;
            if(maxawayU(i) > (relsize*relsize-i) & killpoint == 0)
                killpoint = i
            end
        end
    end
end

figure
plot((maxawayL+maxawayU));

% n actual vs neff
resx = (5*relsize);
resy = (5*relsize);
x = 0:1/resx:Lx;
y = 0:1/resy:Ly;
nodenum = zeros(1,killpoint);
for i = 1:killpoint
    a = zeros(1,(i+maxawayU(i))-(i-maxawayL(i)));
    b = zeros(1,(i+maxawayU(i))-(i-maxawayL(i)));
    psi = zeros(size(x,2),size(y,2));
    psir = zeros(size(x,2),size(y,2));
    psii = zeros(size(x,2),size(y,2));
    for j = (i-maxawayL(i)):(i+maxawayU(i))
        a(j) = randn;
        b(j) = randn;
        for xi = 1:size(x,2)
            for yi = 1:size(y,2)
                psir(xi,yi) = psir(xi,yi) + a(j)*sin(x(xi)*pi*nmtracker(j,1)/Lx)*sin(y(yi)*pi*nmtracker(j,2)/Ly);
                psii(xi,yi) = psii(xi,yi) + b(j)*sin(x(xi)*pi*nmtracker(j,1)/Lx)*sin(y(yi)*pi*nmtracker(j,2)/Ly);
            end
        end
    end
    psi = (psir.*psir + psii.*psii)./(sum(a)*sum(a) + sum(b)*sum(b));
    for xi = 2:(size(x,2)-1)
        for yi = 2:(size(y,2)-1)
            if (psi(xi-1,yi)) < (psi(xi,yi))
                if (psi(xi+1,yi)) < (psi(xi,yi))
                    if (psi(xi,yi-1)) < (psi(xi,yi))
                        if (psi(xi,yi+1)) < (psi(xi,yi))
                            nodenum(i) = nodenum(i) + 1;
                        end
                    end
                end
            end
        end
    end
end
figure
plot(nodenum)
    