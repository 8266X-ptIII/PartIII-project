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
% 

% Written by BGN: 8266X   for part III project
% This builds on the work of Wildeman and utilises his code found in his
% git repo. 


% Dominic's solution implemented with 3 translations 3 rotations in
% cartesian coordinate system
DomSol = fittype('cos(r3)*cos(r2)*(A*((x-t1)*cos(r1)+(y-t2)*sin(r1)).*(besselk(1,sqrt((x-t1).^2+(y-t2).^2)*l1)./(sqrt((x-t1).^2+(y-t2).^2)*besselk(1,R1*l1))))+sin(r2)*y+sin(r3)*x+C', ...
    'independent',{'x','y'}, ...
    'dependent','h', ...
    'problem','R1');


R_in = 24.15; % Radius of the disk, approximate only
pxGrd = 2; % Number of pixels per grid
dist = 35;
l = [38,-1,0.37,0,0.3,0,-R_in,12];
numper = 6; % Number of photos taken per tilt angle.


fold = '2903\Disk\';
folder = strcat('G:/My Drive/Part III Project/',fold,'Set2\');

refname = dir(strcat(folder,'ref*'));
refname = {refname.name};

names = dir(strcat(folder,'*.pgm'));
names = {names.name};
names = names(1:end-length(refname)); % Remove reference images from the photos list.

Iref=imread(strcat(folder,string(refname(1))));
Iref = double(Iref); %[SW]
Iref = flip(Iref,2); % Due to camera oridentation

% get two independent carrier peaks from reference image
[kr, ku] = findorthcarrierpks(Iref, 4*pi/min(size(Iref)), Inf);

% extract carrier signals from reference image and store them for later use
krad = sqrt(sum((kr-ku).^2))/2; %[SW]
fIref = fft2(Iref); %[SW]
cr = getcarrier(fIref, kr, krad); %[SW]
cu = getcarrier(fIref, ku, krad); %[SW]

pxmm = 29.09743116 / (pxGrd*krad) % Calculate the px to mm conversion from lattice spacing

n = 3*numper+1;
h = zeros(size(Iref));
tic
for i=n:n+numper-1 %Itterate over same images at same angle

    name = string(names(i));
    name = char(name);
    angle = str2double(name(1:end-27)) %Extract angle from the file name
    
    Idef=imread(strcat(folder,name));
    
    % convert images to double to prevent rounding errors    
    Idef = double(Idef); %[SW]
    Idef = flip(Idef,2); %Due to camera position

    

    % get displacement field and height profile
    fIdef = fft2(Idef); %[SW]
    [u,v] = fcd_dispfield(fIdef,cr,cu,true);%[SW] % The phase wrap is important
    
    u = u/pxmm;
    v = v/pxmm;    
    %u = medfilt2(u,[5,5]);
    %v = medfilt2(v,[5,5]);
    h = h - mean(h(1:end,end));
    %h = h+invgrad2(-u,-v); % Much Faster    
    h = h+invgrad2(-u,-v)/(pxmm*5.45); %SLOWER BUT PRODUCES MEANINGULLY SCALED RESULTS

end
h=h./numper;

% Rearange for fitting
ys = linspace(1,length(h(:,1)),length(h(:,1)))/pxmm;
xs = linspace(1,length(h(1,:)),length(h(1,:)))/pxmm; 
xss = reshape(repmat(xs,length(ys),1),1,[]).';
yss = repmat(ys(:),length(xs),1);

%h = h(250:end-250,1:end-300);
hss = reshape(h,1,[]).';

[surffit,gof] = fit([xss,yss],hss,DomSol,'Lower',[0,-10,0.2,-0.4,-0.4,-0.5,-R_in], ...
                                        'StartPoint',l, ... % We reuse l as this speeds the process up
                                        'Upper',[1e+4,5,0.5,0.4,0.4,0.3,-R_in], ...
                                        'problem',R_in)

toc

t = num2cell(coeffvalues(surffit));
[A,C,l1,r1,r2,r3,t1,t2] = deal(t{:});


% Prepare data to write to file so can be graphed later
dat = [3;angle;A;l1;gof.rmse;gof.rsquare].';
%writematrix(dat,strcat('G:/My Drive/Part III Project/',fold,'datarotsweakconstraint'),'WriteMode','append');

% Dominic's Solution
mm = linspace(min(xs), max(xs), length(xs));
qq = linspace(min(ys), max(ys), length(ys));
%mm = linspace(min(xs), max(xs), 50);
%qq = linspace(min(ys), max(ys), 50);
[MM, QQ] = meshgrid(mm, qq); 
Z = cos(r2)*cos(r3)*(A*((MM-t1)*cos(r1)+(QQ-t2)*sin(r1)).*(besselk(1,sqrt((MM-t1).^2+(QQ-t2).^2)*l1)./(sqrt((MM-t1).^2+(QQ-t2).^2)*besselk(1,R_in*l1))))+sin(r3)*MM+sin(r2)*QQ+C;


beep
