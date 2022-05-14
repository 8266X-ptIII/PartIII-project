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

% first-order solution
surfexp = fittype('cos(p4)*cos(p3)*(A*exp(((p1*x)+(p2*y))))+sin(p3)*y+sin(p4)*x+C', ...
                    'independent',{'x','y'}, ...
                    'dependent','h');


pxGrd = 2; % Number of pixels per grid

l = [0.0001,0,0.37,0,0,0]; % Initial start point for fitting

fold = ['2903\Razor\Set4\'];
folder = strcat('G:/My Drive/Part III Project/',fold); % folder of photos

refname = dir(strcat(folder,'ref*'));
refname = {refname.name}; % Get reference images saved with the title 'ref'

names = dir(strcat(folder,'*.pgm'));
names = {names.name}; % Get data images
names = names(1:end-length(refname)); % Remove reference images from the photos list.

Iref=imread(strcat(folder,string(refname(3))));
Iref = double(Iref); %[SW]
% 
% % Get two independent carrier peaks from reference image
[kr, ku] = findorthcarrierpks(Iref, 4*pi/min(size(Iref)), Inf); %[SW]
% 
% % Extract carrier signals from reference image and store them for later use
krad = sqrt(sum((kr-ku).^2))/2; %[SW]
fIref = fft2(Iref); %[SW]

cr = getcarrier(fIref, kr, krad); %[SW]
cu = getcarrier(fIref, ku, krad); %[SW]
% 
pxmm = 29.09743116 / (pxGrd*krad) % Calculate the px to mm conversion

%for n = 1:1%:length(names)
n=70;
    name = string(names(n));
    name = char(name);
    angle = str2double(name(1:end-27)) % %27

    Idef=imread(strcat(folder,name));

    %Iref = Idef;
    
    % Get two independent carrier peaks from reference image
    %[kr, ku] = findorthcarrierpks(Iref, 4*pi/min(size(Iref)), Inf); %[SW]

    % Extract carrier signals from reference image and store them for later use
    %krad = sqrt(sum((kr-ku).^2))/2; %[SW]
    %fIref = fft2(Iref); %[SW]
    %cr = getcarrier(fIref, kr, krad); %[SW]
    %cu = getcarrier(fIref, ku, krad); %[SW]

    %pxmm = 29.09743116 / (pxGrd*krad) % Calculate the px to mm conversion
    
    % convert images to double to prevent rounding errors    
    Idef = double(Idef); %[SW]
    
    tic

    % get displacement field and height profile
    fIdef = fft2(Idef); %[SW]
    [u,v] = fcd_dispfield(fIdef,cr,cu,true); %[SW]  % The phase wrap is important
    
    u = u/pxmm;
    v = v/pxmm;    
    %u = medfilt2(u,[5,5]);
    %v = medfilt2(v,[5,5]);
    
    h = invgrad2(-u,-v)/(pxmm*5.45); % Much Faster
    %h = intgrad2(-u,-v,1/pxmm,1/pxmm)/5.45; %SLOWER
    

    % Rearange the data for the reqirements of fit()
    ys = linspace(1,length(h(:,1)),length(h(:,1)))/pxmm;
    xs = linspace(1,length(h(1,:)),length(h(1,:)))/pxmm;    
    xss = reshape(repmat(xs,length(ys),1),1,[]).';
    yss = repmat(ys(:),length(xs),1);


    %figure(4)
    %axis equal
    %quiver(u(1:10:end,1:10:end),v(1:10:end,1:10:end))
    %axis ij
    %axis image

    %figure(5)
    %hold on
    %ax1 = nexttile;
    %mesh(xs,ys,h);
    %set(gca,'DataAspectRatio',[1 1 10])
    %colormap(ax1,hot(8));
    
    % Use one edge to reduce the z offset
    h = h - mean(h(1:end,1));
    hss = reshape(h,1,[]).';
    
    % Perform the surface fitting
     % We reuse l as this speeds the process up
    [surffit,gof] = fit([xss,yss],hss,surfexp','Lower',[0,-25,0.2,-0.1,-2,-1.5],'StartPoint',l,'Upper',[2,10,0.8,0.1,3.0,1])
    l = coeffvalues(surffit);
    conLim = confint(surffit);
    toc

    close all;

    A = l(1);
    l1 = l(3);
    l2 = l(4);
    lc = 1/(l1^2+l2^2)^(1/2)
    errorlc =lc* (abs(conLim(1,3)-l1))/l1;


    %dat = [38;7;angle;lc;errorlc;A].';
    %writematrix(dat,strcat('G:/My Drive/Part III Project/',fold,'datawithnoref'),'WriteMode','append');

    % Rest is to plot a check of the data.

    %mm = linspace(min(xs), max(xs), 50);%length(xs));
    %qq = linspace(min(ys), max(ys), 50);%length(ys));
    mm = linspace(min(xs), max(xs), length(xs));
    qq = linspace(min(ys), max(ys), length(ys));
    [MM, QQ] = meshgrid(mm, qq);
    
    % Exponential Fit
    Z = cos(l(6))*cos(l(5))*(A*exp((MM*l1+QQ*l2)))+sin(l(5))*QQ+sin(l(6))*MM+l(2);
    %figure(6)
    %hold on
    %q=surf(MM, QQ, Z);
    %mesh(xs,ys,h);
      
%end

beep




