clear
randn('state',sum(100*clock));
rand('state',sum(100*clock));

nodenum = zeros(10,1000);
for k = 1:10
    tic
    %defines the size of the box
    wavelength = 1/k; %times 2pi
    boxx = 32;
    boxy = 32;
    xres = 1/k;
    yres = 1/k;
    x = 0:xres:boxx;
    y = 0:yres:boxy;
    %defines how many waves are generated
    waves = 1000;
    %initializes the random variables
    a = zeros(1,waves);
    theta = zeros(1,waves);
    delta = zeros(1,waves);
    %sets the number of times a certian box size is generated
    iter = 100;
    
    for j = 1:iter
        %initializes the space for each iteration
        val = zeros(1,iter);
        psi = zeros(size(x,2),size(y,2));
        psir = zeros(size(x,2),size(y,2));
        psii = zeros(size(x,2),size(y,2));
        for i = 1:waves
            %generates the random components
            a(i) = randn; %iid normal, mean 0, (for real)
            b(i) = randn; %iid normal, mean 0 (for real)
            theta(i) = rand*2*pi; %iid uniform 0 to 2pi
            delta(i) = rand*2*pi; %iid uniform 0 to 2pi
            thetai(i) = rand*2*pi; %iid uniform 0 to 2pi
            deltai(i) = rand*2*pi; %iid uniform 0 to 2pi
        end
        for i = 1:waves
            for xi = 1:size(x,2)
                for yi = 1:size(y,2)
                    %superposition of plane waves
                    psir(xi,yi) = psir(xi,yi) + a(i)/sqrt(sum(a.^2))*cos(x(xi)*cos(theta(i))/wavelength + y(yi)*sin(theta(i))/wavelength + delta(i));
                    psii(xi,yi) = psii(xi,yi) + b(i)*sin(x(xi)*cos(thetai(i))/wavelength + y(yi)*sin(thetai(i))/wavelength + deltai(i));
                end
            end
        end
        psi = psir.*psir;% + psii.^2;
        
%         % start local maximum finder
%         for xi = 2:(size(x,2)-1)
%             for yi = 2:(size(y,2)-1)
%                 if (psi(xi-1,yi)) < (psi(xi,yi))
%                     if (psi(xi+1,yi)) < (psi(xi,yi))
%                         if (psi(xi,yi-1)) < (psi(xi,yi))
%                             if (psi(xi,yi+1)) < (psi(xi,yi))
%                                 nodenum(k,j) = nodenum(k,j) + 1;
%                                 val(j) = val(j) + psi(xi,yi);
%                                 %nodex(nodenum) = x(xi);
%                                 %nodey(nodenum) = y(yi);
%                                 
% %                                 % Full width at half maximum approximation
% %                                 % to find standard deviation of node
% %                                 % First we look in the -x direction
% %                                 for i = 1:(xi-1)
% %                                     if (psi(xi-i,yi))^2 < 0.5*(psi(xi,yi))^2
% %                                         % if a value less than half the max
% %                                         % is found then we interpolate back
% %                                         % linearly to find the FWHM
% %                                         fwhm(nodenum) = 2*abs(abs(x(xi-i+1)-x(xi-i))/2 - nodex(nodenum));
% %                                     else
% %                                         % we make sure that the derivative
% %                                         % stays negative along the vector
% %                                         % of our search
% %                                         if (psi(xi-i-1),yi)^2 > (psi(xi-i),yi)^2
% %                                             i = size(x,2)+1;
% %                                         end
% %                                     end
% %                                 end
% %                                 % If the FWHM has not been found we then
% %                                 % look in the +x direction
% %                                 if fwhm(nodenum) = 0
% %                                     for i = 1:(size(x,2)-xi)
% %                                         if (psi(xi+i,yi))^2 < 0.5*(psi(xi,yi))^2
% %                                             fwhm(nodenum) = 2*abs(abs(x(xi+i-1)-x(xi+i))/2 - nodex(nodenum));
% %                                         else
% %                                             if (psi(xi+i+1),yi)^2 > (psi(xi+i),yi)^2
% %                                                 i = size(x,2)+1;
% %                                             end
% %                                         end
% %                                     end
% %                                 end
% %                                 % If the FWHM has not been found yet then
% %                                 % we look in the -y direction
% %                                 if fwhm(nodenum) = 0
% %                                     for i = 1:(yi-1)
% %                                         if (psi(xi,yi-i))^2 < 0.5*(psi(xi,yi))^2
% %                                             fwhm(nodenum) = 2*abs(abs(y(yi-i+1) - y(yi-i))/2 - nodey(nodenum));
% %                                         else
% %                                             if (psi(xi,yi-i-1))^2 > (psi(xi,yi-i))^2
% %                                                 i = size(y,2)+1;
% %                                             end
% %                                         end
% %                                     end
% %                                 end
% %                                 % If the FWHM has not been found yet then
% %                                 % we look in the +y direction
% %                                 if fwhm(nodenum) = 0
% %                                     for i = 1:(size(y,2)-yi)
% %                                         if (psi(xi,yi+i))^2 < 0.5*(psi(xi,yi))^2
% %                                             fwhm(nodenum) = 2*abs(abs(y(yi+i-1) - y(yi+i))/2 - nodey(nodenum));
% %                                         else
% %                                             if (psi(xi,yi+i+1))^2 > (psi(xi,yi+1))^2
% %                                                 i = size(y,2)+1;
% %                                             end
% %                                         end
% %                                     end
% %                                 end
% %                                 sd(nodenum) = fwhm(nodenum)/2.354;
% %                                 % End of FWHM
%                                 
%                             end
%                         end
%                     end
%                 end
%             end
%         end % end local max finder

        high(j) = max(max(psi)); % extreme value for the set
         
    end
    hightot(k,:) = high;
    %valtot(k,:) = val;
    toc
    
%     % Randomness tests
%     % 1) Values of random points should look gaussian.
%     for i = 1:1000
%         testpx = floor(rand*size(x,2))+1;
%         testpy = floor(rand*size(y,2))+1;
%         val(i) = psi(testpx,testpy);
%     end
%     res = 100;
%     sep = min(val):(1/res)*(max(val)-min(val)):max(val);
%     sepp = zeros(1,res+1);
%     for i = 1:1000
%         for j = 1:100
%             if(val(i) < sep(j))
%                 sepp(j) = sepp(j) + 1;
%                 val(i) = max(val)+1;
%             end
%         end
%     end
%     sepp = sepp*res;
%     figure
%     plot(sep,sepp);
%     % 2) Finds the two point correlation function and makes sure it
%     % approximates a J_0 bessel function.
%     p = zeros(2*size(x,2),2*size(y,2));
%     for i = 0:(2*size(x,2)-1)
%         for j = 0:(2*size(y,2)-1)
%             for m = 0:(size(x,2)-1)
%                 for n = 0:(size(y,2)-1)
%                     p(i+1,j+1) = p(i+1,j+1) + psi(m+1,n+1)*psi(mod(m+i,size(x,2))+1,mod(n+j,size(y,2))+1);
%                 end
%             end
%         end
%     end
%     figure
%     surf(p);
%     % end of randomness tests
    
%     % avg number of local maxima
%     nodetotal(k) = nodenum/iter
end
%nodetotal