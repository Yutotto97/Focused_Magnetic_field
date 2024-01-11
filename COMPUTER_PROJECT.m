%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Yuto
% Show how to plot magnetic fields using the scripts 
% FieldSolenoid & FieldBar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1 solenoid
%clear all; close all; clc;
a=5;        % Solenoid radius
L=3;        % Solenoid length
Br=1;       % Residual field = field in infinite solenoid (mu*n*I)
Theta=0;    % Angle of rotation
RM=[cos(Theta) -sin(Theta);sin(Theta) cos(Theta)];    % Matrice of rotation
[X1,Y1] = meshgrid(-20:5/(3*pi):20);

for k=1:length(X1)
    for kk=1:length(Y1)
        V=[X1(k,kk);Y1(k,kk)];                        % Choosing coordinates
        VR=RM*V;                                    % Rotating the coordinate                        
        BB=FieldSolenoid(a,L,Br,VR(1),VR(2));       % BB(1:2)=(Brho,Bz)
        Brho=BB(1); Bz1(k,kk)=BB(2);
        Bx1(k,kk)=Brho*sign(X1(k,kk)); 
        B1(k,kk)=sqrt(BB(1)^2+BB(2)^2);
    end
end

Ac=[-a ; L/2]; AR=RM*Ac;   %coordinates of the solenoid corners
Bc=[a ; L/2]; BR=RM*Bc;
Cc=[a ; -L/2]; CR=RM*Cc;
Dc=[-a ; -L/2]; DR=RM*Dc;

figure(); contourf(X1,Y1,B1,100,'EdgeColor','none'); colorbar; 

hold on;
plot([AR(1) BR(1)],[AR(2) BR(2)],'-k','Linewidth',2); 
plot([BR(1) CR(1)],[BR(2) CR(2)],'-k','Linewidth',2);
plot([CR(1) DR(1)],[CR(2) DR(2)],'-k','Linewidth',2); 
plot([DR(1) AR(1)],[DR(2) AR(2)],'-k','Linewidth',2);
rl=10;         %proportion of red lines
h = streamline(X1,Y1,Bx1,Bz1,X1(1:rl:end,1:rl:end),Y1(1:rl:end,1:rl:end)); set(h,'Color','red');


%% 3 solenoids
%clear all; clc;
a=5;        % Solenoid radius
L=1;      % Solenoid length
Br=1;       % Residual field  = field in infinite solenoid (mu*n*I)

Theta=pi/6;   % Angle between two most distant magnets
Radius=10;    %radius of partial torus

Center=[Radius;0];
RM1=[cos(Theta/2) -sin(Theta/2);sin(Theta/2) cos(Theta/2)];
RM2=[cos(Theta/2) sin(Theta/2);-sin(Theta/2) cos(Theta/2)];

[X3,Y3] = meshgrid(-20:5/(3*pi):20);

for k=1:length(X3)
    for kk=1:length(Y3)
        V1=[X3(k,kk);Y3(k,kk)];                        % Choosing coordinates
        V1p=V1-Center;
        VR1=RM1*V1p;                                  % Rotating the coordinate
        VR2=RM2*V1p;
        VR1=VR1+Center;
        VR2=VR2+Center;
        BB1=FieldSolenoid(a,L,Br,VR1(1),VR1(2));    
        BB2=FieldSolenoid(a,L,Br,VR2(1),VR2(2));
        BB3=FieldSolenoid(a,L,Br,V1(1),V1(2));
        BB=BB1+BB2+BB3; 
        Brho=BB(1);
        Bx3(k,kk)=BB1(1)*sign(VR1(1))+BB2(1)*sign(VR2(1))+BB3(1)*sign(V1(1));
        Bz3(k,kk)=BB(2); 
        B3(k,kk)=sqrt(BB(1)^2+BB(2)^2);
    end
end

Ac=[-a ; L/2];                 %coordinates of the solenoid corners
Acc=Ac-Center;
AR=RM1*Acc;   
AR=AR+Center;

Bc=[a ; L/2]; Bcc=Bc-Center; BR=RM1*Bcc; BR=BR+Center;

Cc=[a ; -L/2]; Ccc=Cc-Center; CR=RM1*Ccc; CR=CR+Center;

Dc=[-a ; -L/2]; Dcc=Dc-Center; DR=RM1*Dcc; DR=DR+Center;

Ec=Acc; ER=RM2*Ec; ER=ER+Center;

Fc=Bcc; FR=RM2*Fc; FR=FR+Center;

Gc=Ccc; GR=RM2*Gc; GR=GR+Center;

Hc=Dcc; HR=RM2*Hc; HR=HR+Center;



figure(); contourf(X3,Y3,B3,100,'EdgeColor','none'); colorbar; hold on; %colormap('Hot')


plot([AR(1) BR(1)],[AR(2) BR(2)],'-k','Linewidth',2); plot([BR(1) CR(1)],[BR(2) CR(2)],'-k','Linewidth',2);
plot([CR(1) DR(1)],[CR(2) DR(2)],'-k','Linewidth',2); plot([DR(1) AR(1)],[DR(2) AR(2)],'-k','Linewidth',2);

plot([ER(1) FR(1)],[ER(2) FR(2)],'-k','Linewidth',2); plot([FR(1) GR(1)],[FR(2) GR(2)],'-k','Linewidth',2);
plot([GR(1) HR(1)],[GR(2) HR(2)],'-k','Linewidth',2); plot([HR(1) ER(1)],[HR(2) ER(2)],'-k','Linewidth',2);

plot([Ac(1) Bc(1)],[Ac(2) Bc(2)],'-k','Linewidth',2); plot([Bc(1) Cc(1)],[Bc(2) Cc(2)],'-k','Linewidth',2);
plot([Cc(1) Dc(1)],[Cc(2) Dc(2)],'-k','Linewidth',2); plot([Dc(1) Ac(1)],[Dc(2) Ac(2)],'-k','Linewidth',2);

% rl=10;         %proportion of red lines
% h = streamline(X3,Y3,Bx3,Bz3,X3(1:rl:end,1:rl:end),Y3(1:rl:end,1:rl:end)); set(h,'Color','red');

for x=1:length(X3)
    for y=1:length(Y3)
        DB13(x,y)=(B3(x,y)-B1(x,y))/B1(x,y)*100;
    end
end

figure(); contourf(X3,Y3,DB13,100,'EdgeColor','none'); colorbar; hold on; %colormap('Hot')

figure();
plot(X3(38,:),DB13(38,:));

%% 5 solenoids

%clear all; clc;
a=5;        % Solenoid radius
Br=1;       % Residual field  = field in infinite solenoid (mu*n*I)

Theta=pi/6;   % Angle between two most distant magnets
Radius=10;    %radius of partial torus
Nm=5;         %Number of solenoid
L=3/Nm;      % Solenoid length

Center=[Radius;0];
RMs=[];
Thets=[-Theta:2*Theta/(Nm-1):Theta];
for i=1:Nm
    RM=[cos(Thets(i)) -sin(Thets(i));sin(Thets(i)) cos(Thets(i))];
    RMs=[RMs;RM];
end

[X5,Y5] = meshgrid(-20:5/(3*pi):20);

for k=1:length(X5)
    for kk=1:length(Y5)
        BB=0;
        for j=1:2:Nm*2
            V=[X5(k,kk);Y5(k,kk)];                       % Choosing coordinates
            Vc=V-Center;
            VR=RMs(j:j+1,1:2)*Vc;
            VR=VR+Center;
            BBp=FieldSolenoid(a,L,Br,VR(1),VR(2));
            BB=BB+BBp;   
            Bx5(k,kk)=Bx5(k,kk)+BBp(1)*sign(VR(1));
        end 
        Brho=BB(1);
        Bz5(k,kk)=BB(2); 
        B5(k,kk)=sqrt(BB(1)^2+BB(2)^2);
    end
end


figure(); contourf(X5,Y5,B5,100,'EdgeColor','none'); colorbar; hold on; %colormap('Hot')

for q=1:2:Nm*2
    Ac=[-a ; L/2];                 %coordinates of the solenoid corners
    Acc=Ac-Center;
    AR=RMs(q:q+1,1:2)*Acc;   
    AR=AR+Center;
    Bc=[a ; L/2]; Bcc=Bc-Center; BR=RMs(q:q+1,1:2)*Bcc; BR=BR+Center;
    Cc=[a ; -L/2]; Ccc=Cc-Center; CR=RMs(q:q+1,1:2)*Ccc; CR=CR+Center;
    Dc=[-a ; -L/2]; Dcc=Dc-Center; DR=RMs(q:q+1,1:2)*Dcc; DR=DR+Center;
    plot([AR(1) BR(1)],[AR(2) BR(2)],'-k','Linewidth',2); plot([BR(1) CR(1)],[BR(2) CR(2)],'-k','Linewidth',2);
    plot([CR(1) DR(1)],[CR(2) DR(2)],'-k','Linewidth',2); plot([DR(1) AR(1)],[DR(2) AR(2)],'-k','Linewidth',2);
end

% rl=10;         %proportion of red lines
% h = streamline(X5,Y5,Bx5,Bz5,X5(1:rl:end,1:rl:end),Y5(1:rl:end,1:rl:end)); set(h,'Color','red');

for x=1:length(X5)
    for y=1:length(Y5)
        DB15(x,y)=(B5(x,y)-B1(x,y))/B1(x,y)*100;
    end
end

figure(); contourf(X5,Y5,DB15,100,'EdgeColor','none'); colorbar; hold on; %colormap('Hot')

figure();
plot(X5(38,:),DB15(38,:));


%% Prototype initial condition

a=5;        % Solenoid radius
L=1;      % Solenoid length
Br=1;       % Residual field  = field in infinite solenoid (mu*n*I)

Lr=10;  %length of the flexible rod
A_2S=5; %distance between two solenoids at d=0

Alpha=0;

D_2S=[0;A_2S];

[Xp,Yp] = meshgrid(-20:2.5/(3*pi):20);

for k=1:length(Xp)
    for kk=1:length(Yp)
        V1=[Xp(k,kk);Yp(k,kk)];                        % Choosing coordinates
        V1p=V1+D_2S;
        V2p=V1-D_2S;
        BB1=FieldSolenoid(a,L,Br,V1p(1),V1p(2));   
        BB2=FieldSolenoid(a,L,Br,V2p(1),V2p(2));
        BB3=FieldSolenoid(a,L,Br,V1(1),V1(2));
        BB=BB1+BB2+BB3; 
        Brho=BB(1);
        Bxp(k,kk)=BB1(1)*sign(V1p(1))+BB2(1)*sign(V2p(1))+BB3(1)*sign(V1(1));
        Bzp(k,kk)=BB(2); 
        Bp(k,kk)=sqrt(BB(1)^2+BB(2)^2);
    end
end

Ac=[-a ; L/2];                 %coordinates of the solenoid corners
Acc=Ac+D_2S;

Bc=[a ; L/2]; Bcc=Bc+D_2S;

Cc=[a ; -L/2]; Ccc=Cc+D_2S;

Dc=[-a ; -L/2]; Dcc=Dc+D_2S;

Ec=Ac-D_2S; Fc=Bc-D_2S; Gc=Cc-D_2S; Hc=Dc-D_2S;

figure(); 
contourf(Xp,Yp,Bp,100,'EdgeColor','none'); colorbar; hold on;


plot([Acc(1) Bcc(1)],[Acc(2) Bcc(2)],'-k','Linewidth',2); plot([Bcc(1) Ccc(1)],[Bcc(2) Ccc(2)],'-k','Linewidth',2);
plot([Ccc(1) Dcc(1)],[Ccc(2) Dcc(2)],'-k','Linewidth',2); plot([Dcc(1) Acc(1)],[Dcc(2) Acc(2)],'-k','Linewidth',2);

plot([Ec(1) Fc(1)],[Ec(2) Fc(2)],'-k','Linewidth',2); plot([Fc(1) Gc(1)],[Fc(2) Gc(2)],'-k','Linewidth',2);
plot([Gc(1) Hc(1)],[Gc(2) Hc(2)],'-k','Linewidth',2); plot([Hc(1) Ec(1)],[Hc(2) Ec(2)],'-k','Linewidth',2);

plot([Ac(1) Bc(1)],[Ac(2) Bc(2)],'-k','Linewidth',2); plot([Bc(1) Cc(1)],[Bc(2) Cc(2)],'-k','Linewidth',2);
plot([Cc(1) Dc(1)],[Cc(2) Dc(2)],'-k','Linewidth',2); plot([Dc(1) Ac(1)],[Dc(2) Ac(2)],'-k','Linewidth',2);

% rl=10;         %proportion of red lines
% h = streamline(Xp,Yp,Bxp,Bzp,Xp(1:rl:end,1:rl:end),Yp(1:rl:end,1:rl:end)); set(h,'Color','red');


%% prototype: d versus center and angle

%This step permits to determine the link between the center and angle
%The operation is extremely long because the vpa solving keeps 32 digits
%We are going to determine the link, and approximate the results to fasten
%the program

a=5;        % Solenoid radius
L=1;      % Solenoid length
Br=1;       % Residual field  = field in infinite solenoid (mu*n*I)

Lr=10;  %length of the flexible rod
A_2S=5; %distance between two solenoids at d=0

ds=0.5:0.5:5;

for i=1:length(ds)
    syms x
    Rd(i)=vpasolve(sqrt(2*ds(i)/x)+(ds(i)/x)^3/2/(6*sqrt(2)) - Lr/x == 0,x); %wolfram research (approximation for arccos)
    %Rs(i)=vpasolve( acos((x-ds(i))/x) - Lr/x == 0, x);
    Alphad(i)=2*A_2S/Rs(i);
end

Rd
Alphad


%% prototype d=0.5:0.5:5 

a=5;        % Solenoid radius
L=1;      % Solenoid length
Br=1;       % Residual field  = field in infinite solenoid (mu*n*I)

Lr=10;  %length of the flexible rod
A_2S=5; %distance between two solenoids at d=0

ds=0.5:0.5:5;

Rs=[100. 50. 33.33 25. 19.99 16.65 14.25 12.44 11.01 9.85];
%determined previsously and verified with Geogebra (not exceeding 10% of
%error)
Alphas=2*[0.05 0.1 0.15 0.20 0.2501 0.3003 0.3509 0.4019 0.4541 0.5076];


filename1 = 'Prototypemagneticfield_A=5.gif';
filename2 = 'Prototypemagneticflux_A=5.gif';
filename3 = 'Prototype_yield_on_abscissa_A=5.gif';

for i=1:length(ds)  
    fig = figure;
    axis tight manual % this ensures that getframe() returns a consistent size

    Center=[Rs(i);0];
    
    Alpha=Alphas(i);
    
    RM1=[cos(Alpha/2) -sin(Alpha/2);sin(Alpha/2) cos(Alpha/2)];
    RM2=[cos(Alpha/2) sin(Alpha/2);-sin(Alpha/2) cos(Alpha/2)];
    
    for k=1:length(Xp)
        for kk=1:length(Yp)
            V1=[Xp(k,kk);Yp(k,kk)];                     % Choosing coordinates
            V1p=V1-Center;
            VR1=RM1*V1p;                                % Rotating the coordinate
            VR2=RM2*V1p;
            VR1=VR1+Center;
            VR2=VR2+Center;      
            BB1=FieldSolenoid(a,L,Br,VR1(1),VR1(2));    
            BB2=FieldSolenoid(a,L,Br,VR2(1),VR2(2));
            BB3=FieldSolenoid(a,L,Br,V1(1),V1(2));
            BB=BB1+BB2+BB3; 
            Brho=BB(1);
            Bxpd(k,kk)=BB1(1)*sign(VR1(1))+BB2(1)*sign(VR2(1))+BB3(1)*sign(V1(1));
            Bzpd(k,kk)=BB(2); 
            Bpd(k,kk)=sqrt(BB(1)^2+BB(2)^2);
        end
    end
    
    for x=1:length(Xp)
        for y=1:length(Yp)
            DB3p(x,y)=(Bpd(x,y)-Bp(x,y))/Bp(x,y)*100;
        end
    end
    
    plot(Xp(76,:),DB3p(76,:));
    title(['yield on abscissa, d = ' num2str(ds(i)) ', max = ' num2str(round(max(DB3p(76,:))*10)/10) ', min = ' num2str(round(min(DB3p(76,:))*10)/10)])
    ylim([-100;220])
    
    Bmax(i)=max(DB3p(76,:));
    indice_max=find(DB3p(76,:)==Bmax(i));
    Xmax(i)=Xp(76,indice_max);
    Bmin(i)=min(DB3p(76,:));
    indice_min=find(DB3p(76,:)==Bmin(i));
    Xmin(i)=Xp(76,indice_min);
    
    drawnow 
    % Capture the plot as an image 
    frame = getframe(fig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1 
       imwrite(imind,cm,filename3,'gif', 'Loopcount',inf); 
    else 
       imwrite(imind,cm,filename3,'gif','WriteMode','append'); 
    end
    
    Ac=[-a ; L/2];                 %coordinates of the solenoid corners
    Acc=Ac-Center;
    AR=RM1*Acc;   
    AR=AR+Center;
    
    Bc=[a ; L/2]; Bcc=Bc-Center; BR=RM1*Bcc; BR=BR+Center;

    Cc=[a ; -L/2]; Ccc=Cc-Center; CR=RM1*Ccc; CR=CR+Center;

    Dc=[-a ; -L/2]; Dcc=Dc-Center; DR=RM1*Dcc; DR=DR+Center;

    Ec=Acc; ER=RM2*Ec; ER=ER+Center;

    Fc=Bcc; FR=RM2*Fc; FR=FR+Center;

    Gc=Ccc; GR=RM2*Gc; GR=GR+Center;

    Hc=Dcc; HR=RM2*Hc; HR=HR+Center;
    
    contourf(Xp,Yp,Bpd,100,'EdgeColor','none'); colorbar; 
    hold on;

    plot([AR(1) BR(1)],[AR(2) BR(2)],'-k','Linewidth',2); plot([BR(1) CR(1)],[BR(2) CR(2)],'-k','Linewidth',2);
    plot([CR(1) DR(1)],[CR(2) DR(2)],'-k','Linewidth',2); plot([DR(1) AR(1)],[DR(2) AR(2)],'-k','Linewidth',2);

    plot([ER(1) FR(1)],[ER(2) FR(2)],'-k','Linewidth',2); plot([FR(1) GR(1)],[FR(2) GR(2)],'-k','Linewidth',2);
    plot([GR(1) HR(1)],[GR(2) HR(2)],'-k','Linewidth',2); plot([HR(1) ER(1)],[HR(2) ER(2)],'-k','Linewidth',2);

    plot([Ac(1) Bc(1)],[Ac(2) Bc(2)],'-k','Linewidth',2); plot([Bc(1) Cc(1)],[Bc(2) Cc(2)],'-k','Linewidth',2);
    plot([Cc(1) Dc(1)],[Cc(2) Dc(2)],'-k','Linewidth',2); plot([Dc(1) Ac(1)],[Dc(2) Ac(2)],'-k','Linewidth',2);
    title(['3 solenoids, d = ' num2str(ds(i))])
    
    drawnow 
    % Capture the plot as an image 
    frame = getframe(fig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1 
       imwrite(imind,cm,filename1,'gif', 'Loopcount',inf); 
    else 
       imwrite(imind,cm,filename1,'gif','WriteMode','append'); 
    end
    
    contourf(Xp,Yp,Bpd,100,'EdgeColor','none'); colorbar; 
    hold on;

    plot([AR(1) BR(1)],[AR(2) BR(2)],'-k','Linewidth',2); plot([BR(1) CR(1)],[BR(2) CR(2)],'-k','Linewidth',2);
    plot([CR(1) DR(1)],[CR(2) DR(2)],'-k','Linewidth',2); plot([DR(1) AR(1)],[DR(2) AR(2)],'-k','Linewidth',2);

    plot([ER(1) FR(1)],[ER(2) FR(2)],'-k','Linewidth',2); plot([FR(1) GR(1)],[FR(2) GR(2)],'-k','Linewidth',2);
    plot([GR(1) HR(1)],[GR(2) HR(2)],'-k','Linewidth',2); plot([HR(1) ER(1)],[HR(2) ER(2)],'-k','Linewidth',2);

    plot([Ac(1) Bc(1)],[Ac(2) Bc(2)],'-k','Linewidth',2); plot([Bc(1) Cc(1)],[Bc(2) Cc(2)],'-k','Linewidth',2);
    plot([Cc(1) Dc(1)],[Cc(2) Dc(2)],'-k','Linewidth',2); plot([Dc(1) Ac(1)],[Dc(2) Ac(2)],'-k','Linewidth',2);
    title(['3 solenoids, d = ' num2str(ds(i)) ])
    
    rl=10;         
    h = streamline(Xp,Yp,Bxpd,Bzpd,Xp(1:rl:end,1:rl:end),Yp(1:rl:end,1:rl:end)); set(h,'Color','red');
    
    drawnow 
    % Capture the plot as an image 
    frame = getframe(fig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);
    
    if i == 1 
       imwrite(imind,cm,filename2,'gif', 'Loopcount',inf); 
    else 
       imwrite(imind,cm,filename2,'gif','WriteMode','append'); 
    end 
    
    close()
end

figure()
plot(ds,Bmax);
title('minimum yield plotted against d')
figure()
plot(ds,Xmax);
title('position of minimum yield plotted against d')
figure()
plot(Xmax,Bmax);
title('minimum yield plotted against its position')

figure()
plot(ds,Bmin);
title('minimum yield plotted against d')
figure()
plot(ds,Xmin);
title('position of minimum yield plotted against d')
figure()
plot(Xmin,Bmin);
title('minimum yield plotted against its position')





