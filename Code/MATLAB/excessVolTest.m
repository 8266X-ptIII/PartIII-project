% MIT License
% 
% Copyright (c) 2017 Sander Wildeman
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% Lines commented with '[SW]' are attributed to Wildeman see
% https://github.com/swildeman/fcd

% Written by BGN: 8266X  for part III project
% This builds on the work of Wildeman and utilises his code found in his
% git repo. 


gaussian = @(x,s,t) exp(-(x-t).^2/(2*s^2))/(sqrt(2*3.1415)*s);
skewedgaussian = @(x,s,t,alpha) 2*gaussian(x,s,t).*normcdf(alpha*(x-t)/s);


surfgauss = fittype('A*(exp(-(x-t1).^2/(2*s1^2))/(sqrt(2*3.1415)*s1)).*normcdf(w1*(x-t1)/s1).*(exp(-(y-t2).^2/(2*s2^2))/(sqrt(2*3.1415)*s2))','independent',{'x','y'},'dependent','h');

surfjons = fittype('A*exp(-0.5*(w1+asinh((x-t1)/s1)).^2).*exp(-(y-t2).^2/(2*s2^2))./(2*3.1415*s2*s1)+z1','independent',{'x','y'},'dependent','h');



%for set=1:5

fold = ['2903\Razor\Set4\'];
fold = ['1404\Set1\'];
folder = strcat('G:\My Drive\Part III Project\',fold);
%folder = strcat('G:/My Drive/Part III Project/',fold,'Set',string(set),'\');

names = dir(strcat(folder,'*.pgm'));
names = {names.name};
names = names(1:end);%-length(refname)); % Remove reference images from the photos list.
%refname = dir(strcat(folder,'ref*'));
%refname = {refname.name};
%names = refname;


for n=35:35
    close all;

    %Iref = double(imread(strcat(folder,string(refname(n)))));
    name = string(names(n));
    Iref = double(imread(strcat(folder,name)));
    disp(name);    
    name = char(name);

    
    
    % get two independent carrier peaks from reference image
    [kr, ku] = findorthcarrierpks(Iref, 4*pi/min(size(Iref)), Inf); %[SW]
    
    % extract carrier signals from reference image and store them for later use
    krad = sqrt(sum((kr-ku).^2))/2;%[SW]
    fIref = fft2(Iref);%[SW]
    cr = getcarrier(fIref, kr, krad);%[SW]
    cu = getcarrier(fIref, ku, krad);%[SW]
    [rows,cols] = size(Iref);%[SW]
    kxvec = fftshift(kvec(cols));%[SW]
    kyvec = fftshift(kvec(rows));%[SW]
    wr = hann(rows,'periodic');%[SW]
    wc = hann(cols,'periodic');%[SW]
    win2d = wr(:)*wc(:)';%[SW]

    fftIm = fftshift(abs(fft2((Iref-mean(Iref(:))).*win2d)));%[SW]


    pxmm = 29.09743116 / (2*krad);
    dat = [str2double(name(1:end-27)),pxmm,kr,ku];
    %dat = [str2double(name(end-8:end-4)),pxmm,kr,ku];
    %dat = [set,n,pxmm,kr,ku];

    c = cu;
    diffIdxX = 30;
    [val,idxX]=min(abs(kxvec-c.k(1)));
    kxvecShort = kxvec(idxX-fix(diffIdxX):idxX+diffIdxX);
    diffIdxY = 6;
    [val,idxY]=min(abs(kyvec-c.k(2)));
    kyvecShort = kyvec(idxY-diffIdxY:idxY+fix(diffIdxY/1));
    
    fftIm1 = fftIm(idxY-diffIdxY:idxY+fix(diffIdxY/1),idxX-fix(diffIdxX):idxX+diffIdxX);
    
    xss = reshape(repmat(kxvecShort,length(kyvecShort),1),1,[]).';
    yss = repmat(kyvecShort(:),length(kxvecShort),1);
    hss = reshape(fftIm1,1,[]).';
    [surffit,gof] = fit([xss,yss],hss,surfjons,'Lower',[0, 0, 0, c.k(1)-0.5, c.k(2)-0.5, -20, -1e5], ...
                                                'StartPoint',[0, 0.01, 0.01, c.k(1), c.k(2), 0, 0], ...
                                               'Upper',[1e5, 0.05, 0.05, c.k(1)+0.5, c.k(2)+0.5, 20, 1e5])
    coef = coeffvalues(surffit);
    conLim = confint(surffit);
    dat = [dat,coef,abs(conLim(1,6)-coef(6))];

    writematrix(dat,strcat('G:/My Drive/Part III Project/',fold,'jons1'),'WriteMode','append');
    
    %mm = linspace(min(kxvecShort), max(kxvecShort), 50);%length(xs));
    %qq = linspace(min(kyvecShort), max(kyvecShort), 50);%length(ys));
    mm = linspace(min(kxvecShort), max(kxvecShort), length(kxvecShort));
    qq = linspace(min(kyvecShort), max(kyvecShort), length(kyvecShort));
    [MM, QQ] = meshgrid(mm, qq);
      
    A = coef(1);
    s1 = coef(2);
    s2 = coef(3);
    t1 = coef(4);
    t2 = coef(5);
    w1 = coef(6);
    z1 = coef(7);

    %Z = A*(exp(-(MM-t1).^2/(2*s1^2))/(sqrt(2*3.1415)*s1)).*normcdf(alpha*(MM-t1)/s1).*(exp(-(QQ-t2).^2/(2*s2^2))/(sqrt(2*3.1415)*s2));
    Z = A*exp(-0.5*(w1+asinh((MM-t1)/s1)).^2).*exp(-(QQ-t2).^2/(2*s2^2))./(2*3.1415*s2*s1)+z1;
        
   % figure(1)
%     %imagesc(kxvecShort, kyvecShort, fftIm,[0,max(fftIm(:))]);
%     hold on;
%     %q=surf(MM, QQ, Z);
%     
%     %q=surf(kxvecShort, kyvecShort, fftIm1);%(Z-fftIm1)./max(fftIm1(:))*100);
%     %q=surf(kxvecShort, kyvecShort, Z);
%     Ax = gca;
%     %Ax.DataAspectRatio = [diff(get(gca, 'XLim')) diff(get(gca, 'XLim')) 4*diff(get(gca, 'ZLim'))];
%     Ax.ZGrid = 'off';
%     Ax.XGrid = 'off';
%     Ax.YGrid = 'off';
%     Ax.Color = 'none';
%     xlabel('k_x [px^{-1}]');
%     %Ax.XLim = [min(kxvecShort),max(kxvecShort)];
%     %Ax.YLim = [min(kyvecShort),max(kyvecShort)];
% 
%     ylabel('k_y [px^{-1}]');
%     %zlabel('Difference [% max height]');
%     zlabel('FFT');
%     Ax.FontName = 'arial';
%     Ax.FontSize = 10;
%     Ax.LineWidth = 1.25;
%     %Ax.ColorScale = 'log';
%     text.color = 'black';
%     a=colorbar();
%     %a.Ticks = [0.05,0.1,0.2,0.5,1,2,3];
%     %a.Ticks = [0.01,0.03,0.06,0.10,0.30,1.00];
%     %caxis([0.05 3]);
%     %caxis([0.01 1.5]);
%     a.Label.String = 'Difference [% max height]';
%     a.Label.String = 'FFT';
%     %a.Label.String = 'z (approximate) [mm]';
%     %view(0,0);
%     %figure_handle = gcf;
%     %cbar_handle = findobj(figure_handle,'tag','Colorbar')
%     %cbar_handle.YAxisLocation = 'right';
%     set(gcf,'units','inches','position',[1,1,10,5])

    dat = [kxvecShort.',max(fftIm1,[],1).',max(Z,[],1).'];

    writematrix(dat,strcat('G:/My Drive/Part III Project/skew plots'),'WriteMode','append');

    figure(10)
    hold on;
    a=area(kxvecShort, max(fftIm1,[],1));
    q=area(kxvecShort, max(Z,[],1));
    
    
    q.FaceAlpha = 0.5;
    q.FaceColor = 'r';
    a.FaceAlpha = 1;
    a.FaceColor = 'g';
    a.LineStyle = 'none';
    q.LineStyle = 'none';
    xlabel('k_x [px^{-1}]');
    ylabel('FFT');
    Ax = gca;
    Ax.FontName = 'arial';
    Ax.FontSize = 10;
    Ax.LineWidth = 1.25;
    set(gcf,'units','inches','position',[1,1,6.75,5])

    %pause(3)

%axis image
end
%end
beep



    