clear
randn('state',sum(100*clock));

load longrunbill.mat;

% relsize = 1500;
% Lx = 1.0;
% Ly = 1.61803399;
% nmax = relsize;
% mmax = floor(relsize*sqrt(Ly));
% 
% E = zeros(nmax,mmax);
% k = 0;
% for n = 1:nmax
%     for m = 1:mmax
%         E(n,m) = (n*n/Lx + m*m/Ly);
%         k = k + 1;
%         nmtracker(k,:) = [n,m];
%     end
% end
% 
% % Sorts the energies, but keeps track of corresponding n,m
% Elist = [];
% for i = 1:nmax
%     Elist = [Elist,E(i,:)];
% end
% 
% for i = 1:(nmax*mmax)
%     currhigh = 0;
%     Ehigh = 0;
%     for j = 1:(nmax*mmax - i + 1)
%         if (Elist(j) >= Ehigh)
%             currhigh = j;
%             Ehigh = Elist(j);
%         end
%     end
%     if currhigh > 0
%         temp = Elist(nmax*mmax - i + 1);
%         Elist(nmax*mmax - i + 1) = Elist(currhigh);
%         Elist(currhigh) = temp;
%         temp2 = nmtracker(nmax*mmax - i + 1,:);
%         nmtracker(nmax*mmax - i + 1,:) = nmtracker(currhigh,:);
%         nmtracker(currhigh,:) = temp2;
%     end
% end
% 
% for i = 1:(nmax*mmax)
%     deltaE(i) = sqrt(2*Elist(i))/((Lx+Ly)/2);
% end
% 
% % plot of E with an error bar of size delta E
% errorbar(Elist,deltaE,'o')
% 
% % neff finder
% killpoint = 0;
% for i = 1:(nmax*mmax)
%     for j = 1:i
%         if((Elist(i) - Elist(j)) < deltaE(i))
%             maxawayL(i) = i-j;
%             break;
%         end
%     end
%     for j = i:(nmax*mmax)
%         if((Elist(j) - Elist(i)) < deltaE(i))
%             maxawayU(i) = j-i;
%             if(maxawayU(i) > (relsize*relsize-i) & killpoint == 0)
%                 killpoint = i
%             end
%         end
%     end
% end
% 
% figure
% plot((maxawayL+maxawayU));

% n actual vs neff
iter = 1000;

for p = 1:10
    resx = 1000;
    resy = 1000;
    x = 0:1/resx:Lx;
    y = 0:1/resy:Ly;
    tic
    i = 2^p;
    nodenum = zeros(1,iter);
    val = zeros(1,iter);
    
    for k = 1:iter

    a = zeros(1,(i+maxawayU(i))-(i-maxawayL(i)));
    b = zeros(1,(i+maxawayU(i))-(i-maxawayL(i)));
    psi = zeros(size(x,2),size(y,2));
    psir = zeros(size(x,2),size(y,2));
    %psii = zeros(size(x,2),size(y,2));
    for j = (i-maxawayL(i)):(i+maxawayU(i))
        a(j) = randn;
        %b(j) = randn;
        for xi = 1:size(x,2)
            for yi = 1:size(y,2)
                psir(xi,yi) = psir(xi,yi) + a(j)*sin(x(xi)*pi*nmtracker(j,1)/Lx)*sin(y(yi)*pi*nmtracker(j,2)/Ly);
                %psii(xi,yi) = psii(xi,yi) + b(j)*sin(x(xi)*pi*nmtracker(j,1)/Lx)*sin(y(yi)*pi*nmtracker(j,2)/Ly);
            end
        end
    end
    psi = (psir.*psir/(sqrt(sum(a.^2))));% + psii.*psii/(sqrt(sum(b.^2))))/2;
    for xi = 2:(size(x,2)-1)
        for yi = 2:(size(y,2)-1)
            if (psi(xi-1,yi)) < (psi(xi,yi))
                if (psi(xi+1,yi)) < (psi(xi,yi))
                    if (psi(xi,yi-1)) < (psi(xi,yi))
                        if (psi(xi,yi+1)) < (psi(xi,yi))
                            nodenum(k) = nodenum(k) + 1;
                            val(k) = val(k) + psi(xi,yi);
                        end
                    end
                end
            end
        end
    end
    maxval(k) = max(max(psi));
    end
    maxtotal(p,:) = maxval;
    nodetotal(p,:) = nodenum;
    nodeavg(p) = sum(nodenum)/1000.0;
    valtotal(p,:) = val;
    valavg(p) = sum(val);
    toc
end

figure
plot(nodenum)
    