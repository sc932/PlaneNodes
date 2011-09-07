clear
randn('state',sum(100*clock));
rand('state',sum(100*clock));

for k = 1:25
    tic
    %defines the size of the box
    wavelength = 1.0;
    boxx = 32;
    boxy = 32;
    xres = 0.5*wavelength;
    yres = 0.5*wavelength;
    x = 0:xres:boxx;
    y = 0:yres:boxy;
    %defines how many waves are generated
    waves = 1000;
    %initializes the random variables
    a = zeros(1,waves);
    theta = zeros(1,waves);
    delta = zeros(1,waves);
    %sets the number of times a certian box size is generated
    iter = 1000;
    nodenum = 0;

    for j = 1:iter
        %initializes the space for each iteration
        psi = zeros(size(x,2),size(y,2));
        for i = 1:waves
            %generates the random components
            a(i) = randn; %iid normal, mean 0, 
            theta(i) = rand*2*pi; %iid uniform 0 to 2pi
            delta(i) = rand*2*pi; %iid uniform 0 to 2pi
            for xi = 1:size(x,2)
                for yi = 1:size(y,2)
                    %superposition of plane waves
                    psi(xi,yi) = psi(xi,yi) + a(i)*sin(x(xi)*cos(theta(i))/wavelength + y(yi)*sin(theta(i))/wavelength + delta(i));
                end
            end
        end
        psi = psi.*psi;
%         for xi = 2:(size(x,2)-1)
%             for yi = 2:(size(y,2)-1)
%                 if (psi(xi-1,yi))^2 < (psi(xi,yi))^2
%                     if (psi(xi+1,yi))^2 < (psi(xi,yi))^2
%                         if (psi(xi,yi-1))^2 < (psi(xi,yi))^2
%                             if (psi(xi,yi+1))^2 < (psi(xi,yi))^2
%                                 nodenum = nodenum + 1;
%                                 yi = yi+1;
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
%        end
        high(j) = max(max(abs(psi))); % extreme value for the set
    end   
    hightotal(k,:) = high
    toc
end