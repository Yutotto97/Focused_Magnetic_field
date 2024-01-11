clear all
clc

I=3;                    % Intensity (A)
nuo=1.26e-6;             % Magnetic constant
d = 5;                  % Radius of the Loop.. (cm)
m = 1;                 % Direction of current through the Loop... 1(ANTI-CLOCK), -1(CLOCK)
N = 30;                 % Number sections/elements in the conductor...
T = 10;                  % Number of Turns in the Solenoid..
g = 0.4;               % Gap between Turns..
dl = 2*pi*d/N;          % Length of each element..
% XYZ Coordinates/Location of each element from the origin(0,0,0), Center of the Loop is taken as origin..
dtht = 360/N;                           % Loop elements divided w.r.to degrees..
tht = (0+dtht/2): dtht : (360-dtht/2);  % Angle of each element w.r.to origin..
xC =  d.*cosd(tht).*ones(1,N);      % X coordinate...
yC =  d.*sind(tht).*ones(1,N);      % Y coordinate...
% zC will change for each turn, -Z to +Z, and will be varied during iteration...
zC = ones(1,N);
h = (d-g*(T-1))/d : g/d : 1;
% Length(Projection) & Direction of each current element in Vector form..
Lx = dl.*cosd(90+tht);        % Length of each element on X axis..
Ly = dl.*sind(90+tht);        % Length of each element on Y axis..
Lz = zeros(1,N);              % Length of each element is zero on Z axis..
% Points/Locations in space (here XZ plane) where B is to be computed..
NP = 125;               % Detector points..
xPmax = 3*d;            % Dimensions of detector space.., arbitrary..
zPmax = 1.5*(T*g);
xP = linspace(-xPmax,xPmax,NP);        % Divide space with NP points..
zP = linspace(-zPmax,zPmax,NP);
[xxP zzP] = meshgrid(xP,zP);            % Creating the Mesh..
% Initialize B..
Bx = zeros(NP,NP);
By = zeros(NP,NP);
Bz = zeros(NP,NP);
% Computation of Magnetic Field (B) using Superposition principle..
% Compute B at each detector points due to each small cond elements & integrate them..
for p = 1:T
    for q = 1:N
        rx = xxP - h(p)*xC(q);               % Displacement Vector along X direction from Conductor..
        ry =  h(p)*yC(q);                   % Displacement Vector along Y direction from Conductor..
        rz = zzP - zC(q);                    % No detector points on Y direction..
    
        r = sqrt(rx.^2+ry.^2+rz.^2);    % Displacement Magnitude for an element on the conductor..
    
        r3 = r.^3;
    
        %        dl X r
        % B = K  ------
        %         r^3
        %                  |  i    j    k   |         |  i    j    k  |
        %   dl X r  =      | Lx    Ly   Lz  |     =   | Lx   Ly    0  |  =   i(Ly.rz-0) - j(Lx.rz-0) + k(Lx.ry - Ly.rx)
        %                  | rx    ry   rz  |         | rx   ry    rz |
    
        % Hence, Bx = Ly.rz/r^3, By = Lx.rz/r^3, Bz = (Lx.ry-Ly.rx)/r^3
        Bx = Bx + nuo*I*m*Ly(q).*rz./r3;      % m - direction of current element..
        By = By - nuo*I*m*Lx(q).*rz./r3;
        Bz = Bz + nuo*I*m*Lx(q).*ry./r3 - nuo*I*m*Ly(q).*rx./r3;
    end
end
B = sqrt(Bx.^2 + By.^2 + Bz.^2);        % Magnitude of B..
%B = B/max(max(B));                      % Normalizing...

% Plotting...
figure();
pcolor(xxP,zzP,B);
colormap(jet);
shading interp;
axis equal;
axis([-xPmax xPmax -zPmax zPmax]);
xlabel('<-- x -->');ylabel('<-- z -->');
title('Magnetic Field Distibution');
colorbar;

% figure();
% surf(xxP,zzP,B,'FaceColor','interp',...
%     'EdgeColor','none',...
%     'FaceLighting','phong');
% daspect([1 1 1]);
% axis tight;
% view(-10,30);
% camlight right;
% colormap(jet);
% grid off;
% axis off;
% colorbar;
% title('Magnetic Field Distibution');
%  
figure();
quiver(xxP,zzP,Bx,Bz);
colormap(lines);
axis([-2*d 2*d -T*g T*g]);
title('Magnetic Field Distibution');
xlabel('<-- x -->');
ylabel('<-- z -->');
zoom on;