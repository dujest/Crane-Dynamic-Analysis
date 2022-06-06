clc
clear
tStart=tic;

%% MODAL BEHAVIOUR
%% Geometrical and material properties of the beams
%%
rho = 7800;      %kg/m^3 - density
E = 210e9;       %Pa - modulus of elasticity (Young's modulus)
v = 0.29;        %Poisson's ratio
G = E/(2*(1+v)); %Pa - shear modlus of elasticity

A1 = 0.0025;     %m^2 - cross-sectional area of beam 1 (L = 1 m)
A2 = 0.0025;     %m^2 - cross-sectional area of beam 2 (L = 2.24 m)
A3 = 0.005;      %m^2 - cross-sectional area of beam 3 (L = 2 m)
A = [A1, A2, A3];

I =zeros(1,length(A)); %m^4 - second moment of area about the axis y and z 
J =zeros(1,length(A)); %m^4 - torsional constant of the cross-section
for i=1:length(A)
    [I(i), J(i)] = BeamProperties(A(i));
end

r=sqrt(I./A);    %radius of gyration of the cross-section

%% Crane jib modelling
%%
Nnode = 21;               %number of nodes

node_list=zeros(Nnode,3); %list of coordinates of each node

node_list(1,:)=[0 0 0];
node_list(2,:)=[1 0 0];
node_list(3,:)=[cosd(60) 0 sind(60)];

yd = 2; %distance between nodes on y axis
for i=3:3:(Nnode-3)
    node_list(i+1,:)=[0 yd 0];
    node_list(i+2,:)=[1 yd 0];
    node_list(i+3,:)=[cosd(60) yd sind(60)];
    
    yd = yd + 2;
end

for i=1:Nnode
    plot3(node_list(i,1),node_list(i,2),node_list(i,3),'.y','markersize',40); %3D plot of the nodes
    text(node_list(i,1),node_list(i,2),node_list(i,3),num2str(i),'horizontalalignment','center','verticalalignment','middle')
    hold on
end

Ncell = 6; %number of cells

beam1_list = zeros((Ncell+1)*3,2);   %list of beams 1 (L = 1 m)
beam1_list(1:3,:) = [1 2; 1 3; 2 3]; %nodes of beam 1 for cell 1

d = 3; %difference between the nodes of each cell
for i=3:3:(Nnode-3)
    beam1_list(i+1,:) = [d+1 d+2];
    beam1_list(i+2,:) = [d+1 d+3];
    beam1_list(i+3,:) = [d+2 d+3];

    d = d + 3;
end

beam2_list = zeros(Ncell*3,2);       %list of beams 2 (L = 2.24 m)
beam2_list(1:3,:) = [2 4; 3 4; 3 5]; %nodes of beam 2 for cell 1

d = 3;
for i=3:3:(Nnode-6)
    beam2_list(i+1,:) = [d+2 d+4];
    beam2_list(i+2,:) = [d+3 d+4];
    beam2_list(i+3,:) = [d+3 d+5];

    d = d + 3; 
end

beam3_list = zeros(Ncell*3,2);       %list of beams 3 (L = 2 m)
beam3_list(1:3,:) = [1 4; 2 5; 3 6]; %nodes of beam 3 for cell 1

d = 3;
for i=3:3:(Nnode-6)
    beam3_list(i+1,:) = [d+1 d+4];
    beam3_list(i+2,:) = [d+2 d+5];
    beam3_list(i+3,:) = [d+3 d+6];

    d = d + 3;
end

n_beam_type = 3; %number of beam types
beam_list = cell(1,n_beam_type);
beam_list{1} = beam1_list;
beam_list{2} = beam2_list;
beam_list{3} = beam3_list;

Nbeam = [length(beam1_list), length(beam2_list), length(beam3_list)]; %number of beams per each type

Color = ['b' 'b' 'g'];
for j=1:n_beam_type
    for i=1:Nbeam(j)
        x = [node_list(beam_list{j}(i,1),1), node_list(beam_list{j}(i,2),1)];
        y = [node_list(beam_list{j}(i,1),2), node_list(beam_list{j}(i,2),2)];
        z = [node_list(beam_list{j}(i,1),3), node_list(beam_list{j}(i,2),3)];
        plot3(x,y,z, 'Color', Color(j))
        hold on
    end
end
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

Nbeams = sum(Nbeam);      %number of beams

L = zeros(1,n_beam_type); %length of each beam type
for i=1:n_beam_type
    L(i) = BeamElementLength(node_list(beam_list{i}(1,1),1), node_list(beam_list{i}(1,1),2), node_list(beam_list{i}(1,1),3), node_list(beam_list{i}(1,2),1), node_list(beam_list{i}(1,2),2), node_list(beam_list{i}(1,2),3));
end

%% Crane jib discretization
%%
n_elem_b1 = 3; %number of elements per beam 1
n_elem_b2 = 3; %number of elements per beam 2
n_elem_b3 = 3; %number of elements per beam 3
n_elem_b = [n_elem_b1, n_elem_b2, n_elem_b3];

Nelements = sum(Nbeam.*n_elem_b);           %number of elements

Nnode_new = Nnode+sum(Nbeam.*(n_elem_b-1)); %new number of nodes
node_list_new = zeros(Nnode_new,3);         %new list of nodes
node_list_new(1:length(node_list),:) = node_list;
element_list = cell(1,3);                   %lists of elements

p = Nnode+1;
for i=1:n_beam_type
    ep = 1;
    element_list{i} = zeros(Nbeam(i)*(n_elem_b(i)),2);
    for j=1:Nbeam(i)
        dx0 =(node_list(beam_list{i}(j,2),1)-node_list(beam_list{i}(j,1),1))/n_elem_b(i);
        dy0 =(node_list(beam_list{i}(j,2),2)-node_list(beam_list{i}(j,1),2))/n_elem_b(i);
        dz0 =(node_list(beam_list{i}(j,2),3)-node_list(beam_list{i}(j,1),3))/n_elem_b(i);
        dx = dx0;
        dy = dy0;
        dz = dz0;
        for k=1:(n_elem_b(i)-1)
            node_list_new(p,:) = [node_list(beam_list{i}(j,1),1)+dx, node_list(beam_list{i}(j,1),2)+dy, node_list(beam_list{i}(j,1),3)+dz];
            if k == 1
                element_list{i}(ep,:) = [beam_list{i}(j,1), p];
                element_list{i}(ep+1,:) = [p, p+1];
                ep = ep+1;
            elseif k == (n_elem_b(i)-1)
                element_list{i}(ep,:) = [p, beam_list{i}(j,2)];
            else
                element_list{i}(ep,:) = [p, p+1];
            end
            dx = dx+dx0;
            dy = dy+dy0;
            dz = dz+dz0;
            p = p+1;
            ep = ep+1;
        end
    end
end

%% Plotting of the mesh
%%
figure
Color = ['b' 'b' 'g'];
for j=1:n_beam_type
    for i=1:Nbeam(j)
        x = [node_list(beam_list{j}(i,1),1), node_list(beam_list{j}(i,2),1)];
        y = [node_list(beam_list{j}(i,1),2), node_list(beam_list{j}(i,2),2)];
        z = [node_list(beam_list{j}(i,1),3), node_list(beam_list{j}(i,2),3)];
        plot3(x,y,z, 'Color', Color(j))
        hold on
    end
end
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
for i=1:Nnode_new
    plot3(node_list_new(i,1),node_list_new(i,2),node_list_new(i,3),'.r','markersize',10); %3D plot of the elements
    hold on
end

%% Assembly of the system matrices
%%
Nelement = [length(element_list{1}), length(element_list{2}), length(element_list{3})]; %number of elements per each beam type

kel = cell(max(Nelement),3); %cell with stiffness matrix of each element 

Ndof = 6; %number of degrees of freedom per node

Nnel = 2; %number of nodes per element

Lel = zeros(n_beam_type,1);

for i=1:n_beam_type
    for j=1:Nelement(i)
        kel{j,i}=zeros(Ndof*Nnel);
        
        Lel(i)=BeamElementLength(node_list_new(element_list{i}(j,1),1), node_list_new(element_list{i}(j,1),2), node_list_new(element_list{i}(j,1),3), node_list_new(element_list{i}(j,2),1), node_list_new(element_list{i}(j,2),2), node_list_new(element_list{i}(j,2),3));
        
        kel{j,i}=BeamElementStiffness(E,G,A(i),I(i),I(i),J(i),node_list_new(element_list{i}(j,1),1), node_list_new(element_list{i}(j,1),2), node_list_new(element_list{i}(j,1),3), node_list_new(element_list{i}(j,2),1), node_list_new(element_list{i}(j,2),2), node_list_new(element_list{i}(j,2),3), Lel(i));
    end
end

K=zeros(Nnode_new*Ndof);

for i=1:n_beam_type
    for j=1:Nelement(i)
        K=CraneJibAssemble(kel{j,i},K,element_list{i}(j,1),element_list{i}(j,2));
    end
end

mel = cell(max(Nelement),3); %cell with mass matrix of each element

for i=1:n_beam_type
    for j=1:Nelement(i)
        mel{j,i}=zeros(Ndof*Nnel);
        
        Lel(i)=BeamElementLength(node_list_new(element_list{i}(j,1),1), node_list_new(element_list{i}(j,1),2), node_list_new(element_list{i}(j,1),3), node_list_new(element_list{i}(j,2),1), node_list_new(element_list{i}(j,2),2), node_list_new(element_list{i}(j,2),3));
        
        mel{j,i}=BeamElementMassMatrix(rho,A(i),node_list_new(element_list{i}(j,1),1), node_list_new(element_list{i}(j,1),2), node_list_new(element_list{i}(j,1),3), node_list_new(element_list{i}(j,2),1), node_list_new(element_list{i}(j,2),2), node_list_new(element_list{i}(j,2),3), Lel(i), r(i));
    end
end

M = zeros(Nnode_new*Ndof); %global mass matrix

for i=1:n_beam_type
    for j=1:Nelement(i)
        M=CraneJibAssemble(mel{j,i},M,element_list{i}(j,1),element_list{i}(j,2)); 
    end
end

%% Boundary conditions
%%
BC = 'clamped';               %'free' or 'clamped'
switch BC
    case 'free'         
        Fixed = [];
    case 'clamped'
        Fixed = [1 2 3];      %clamped nodes
end

fixedDOF = [];                %constrained degrees of freedom
for i=1:length(Fixed)
    fixedDOF = [fixedDOF (Fixed(i)*Ndof-Ndof+1):1:(Fixed(i)*Ndof)];
end

% fixedDOF = 1:18;
NfixedDOF = length(fixedDOF); %number of constrained dofs

K(fixedDOF, :) = [];
K(:, fixedDOF) = [];
M(fixedDOF, :) = [];
M(:, fixedDOF) = [];

%% Eigenvalues problem
%%
[V, D] = eigs(K,M,10,'smallestabs');

w = sqrt(diag(D));   %rad/s - first 10 angular eigenfrequencies
f = w/(2*pi);        %Hz - first 10 eigenfrequencies

Neigenfrequency = 7; %number of eigenfrequencies
tEigenfrequencies=toc(tStart);

%% Plotting of the eigenmodes
%%
deformed_node_list = node_list_new; %coordinates of the deformed nodes
scale = 15;                         %scale of the deformed structure
Neigenmode = 3;                     %the first three eigenmodes

for n = 1:Neigenmode
    for i=4:Nnode_new
        deformed_node_list(i,:) = node_list_new(i,:)+scale*V((1+Ndof*(i-3)-Ndof):(Ndof*(i-3)-3),n)';
    end
    
    figure
    Color = ['b' 'b' 'g'];
    for j=1:n_beam_type
        %plotting of the undeformed structure
        for i=1:Nbeam(j)
            x = [node_list(beam_list{j}(i,1),1), node_list(beam_list{j}(i,2),1)];
            y = [node_list(beam_list{j}(i,1),2), node_list(beam_list{j}(i,2),2)];
            z = [node_list(beam_list{j}(i,1),3), node_list(beam_list{j}(i,2),3)];
            plot3(x,y,z, ':k', 'LineWidth', 1)
            hold on
        end
        
        %plotting of the deformed structure
        for i=1:Nelement(j)
            x = [deformed_node_list(element_list{j}(i,1),1), deformed_node_list(element_list{j}(i,2),1)];
            y = [deformed_node_list(element_list{j}(i,1),2), deformed_node_list(element_list{j}(i,2),2)];
            z = [deformed_node_list(element_list{j}(i,1),3), deformed_node_list(element_list{j}(i,2),3)];
            plot3(x,y,z, '-', 'LineWidth', 1, 'Color', Color(j))
            hold on
        end
    end
    axis equal
    title(['Eigenmode ',num2str(n),' (',num2str(f(n)),' Hz)'])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
end

%% Damping Ratios
%%
damp_ratio = zeros(Neigenfrequency,1);   %number of damping ratios
damp_ratio(1) = 0.006;
damp_ratio(2) = 0.008;

a = (2*damp_ratio(1)*w(1)-2*damp_ratio(2)*w(2))/((w(1))^2-(w(2))^2); %stiffness damping constant
b = 2*damp_ratio(1)*w(1)-a*((w(1))^2);                               %mass damping constant

C = a*K+b*M;                             %damping matrix

for i=3:Neigenfrequency
    damp_ratio(i) = 0.5*(a*w(i)+b/w(i)); %first 7 damping ratios
end

%% DYNAMIC RESPONSE TO AN EXTERNAL LOAD
%% Superpostion Methods
%%
dt = 0.01;         %s - time step
time = 1;         %s - total time
t = 0:dt:time;     %time interval
Nstep = length(t); %number of steps

T = 0.7;           %s - harmonic excitation period
%T = 0.2;
F0 = 400;          %N - harmonic exication force

Ndofs = length(K); %number of free degrees of freedom
F = zeros(Ndofs,Nstep);

nHF = [19 20];     %nodes with harmonic force in x direction

for i=1:length(nHF)
    F(nHF(i)*Ndof-Ndof+1-18,:) = F0*sin((2*pi/T)*t); %applied harmonic forces
end

Gmass = zeros(Neigenfrequency,1);     %generalized masses
wD    = zeros(Neigenfrequency,1);     %damped angular eigenfrequencies
phi   = zeros(Neigenfrequency,Nstep); %modal participation factor
h     = zeros(Neigenfrequency,Nstep); 
eta   = zeros(Neigenfrequency,Nstep); 

for i=1:Neigenfrequency
    Gmass(i,:) = V(:,i)'*M*V(:,i);
    wD(i,:) = w(i)*sqrt(1-(damp_ratio(i))^2);
    
    for j=1:Nstep
        phi(i,j) = V(:,i)'*F(:,j)/Gmass(i);    
        
        h(i,j) = (dt/wD(i))*(exp(-damp_ratio(i)*w(i)*t(j))*sin(wD(i)*t(j)));
    end
    
    convolution = conv(phi(i,:),h(i,:));
    eta(i,:) = convolution(1:j);
end

%% Modal displacement Method
%%
q = zeros(Ndofs,Nstep);          %respones of the jib

for j=1:Nstep 
    for i=1:Neigenfrequency
        q(:,j) = q(:,j)+eta(i,j)*V(:,i);
    end    
end

%% Modal Acceleration Method
%%
third_term = zeros(Ndofs,Nstep); %third term of the response equation
q2 = zeros(Ndofs,Nstep);         %response of the jib

for j=1:Nstep 
    for i=1:Neigenfrequency
        third_term(:,j) = third_term(:,j)+(phi(i,j)/(w(i)*w(i)))*V(:,i);
        
        q2(:,j) = (q(:,j)+eta(i,j)*V(:,i))+((K^-1)*F(:,j)-third_term(:,j));
    end    
end

%% Newmark Integration Scheme
%%
q0    = zeros(Ndofs,1);         %displacement at 0
q0I   = zeros(Ndofs,1);         %velocity at 0
q0II  = M^-1*F(:,1)-C*q0I-K*q0; %acceleration at 0

q3(1,:)  = q0;
qI(1,:)  = q0I;
qII(1,:) = q0II;

q_star(1,:)   = q0;
qI_star(1,:)  = q0I;
qII_star(1,:) = q0II;

h2 = 0.01;      %s - time increment
t2 = 0:h2:time; %s - total time
gamma = 0.5;
beta  = 0.25;

S = M+gamma*h2*C+h2^2*beta*K;

for i=1:length(t2)-1
    %Prediction
    qI_star(i+1,:) = qI(i,:)+(1-gamma)*h2*qII(i,:);
    q_star(i+1,:) = q3(i,:)+h2*qI(i,:)+(0.5-beta)*h2^2*qII(i,:);
    
    %Acceleration calculation
    qII(i+1,:) = S^-1*(F(:,i+1)-C*qI_star(i+1,:)'-K*q_star(i+1,:)');
    
    %Correction
    qI(i+1,:) = qI_star(i+1,:)+h2*gamma*qII(i+1,:);
    q3(i+1,:) = q_star(i+1,:)+h2^2*beta*qII(i+1,:);    
end

%% Plotting of the dynamic responses
%%
figure
for i=1:length(nHF)
     plot(t,F(nHF(i)*Ndof-Ndof+1-18,:)); %applied harmonic forces
     hold on
end
title('Applied Forces')
xlabel('Time [s]')
ylabel('Force [N]')

%Displacement of the node A along x direction
nodeA = 20; %number of the node A in the node list
figure
subplot(4,1,1)
plot(t,q(nodeA*Ndof-Ndof+1-18,:),'-b')
title('Dynamic Response - X Direction')
xlabel(['Time [s] (Time increment = ', num2str(dt),' s)'])
ylabel('Displacement [m]')
legend('Modal Displacement Method')

subplot(4,1,2)
plot(t,q2(nodeA*Ndof-Ndof+1-18,:),'-g')
xlabel(['Time [s] (Time increment = ', num2str(dt),' s)'])
ylabel('Displacement [m]')
legend('Modal Acceleration Method')

subplot(4,1,3)
plot(t2,q3(:,nodeA*Ndof-Ndof+1-18)','-r')
xlabel(['Time [s] (Time increment = ', num2str(h2),' s)'])
ylabel('Displacement [m]')
legend('Newmark Integration Scheme')

subplot(4,1,4)
plot(t2,q(nodeA*Ndof-Ndof+1-18,:),'b',t,q2(nodeA*Ndof-Ndof+1-18,:),'g',t,q3(:,nodeA*Ndof-Ndof+1-18)','r')
title('X Direction Comparisson')
xlabel(['Time [s] (Time increment = ', num2str(h2),' s)'])
ylabel('Displacement [m]')

%Displacement of the node A along y direction
figure
subplot(4,1,1)
plot(t,q(nodeA*Ndof-Ndof+2-18,:),'-b')
title('Dynamic Response - Y Direction')
xlabel(['Time [s] (Time increment = ', num2str(dt),' s)'])
ylabel('Displacement [m]')
legend('Modal Displacement Method')

subplot(4,1,2)
plot(t,q2(nodeA*Ndof-Ndof+2-18,:),'-g')
xlabel(['Time [s] (Time increment = ', num2str(dt),' s)'])
ylabel('Displacement [m]')
legend('Modal Acceleration Method')

subplot(4,1,3)
plot(t2,q3(:,nodeA*Ndof-Ndof+2-18)','-r')
xlabel(['Time [s] (Time increment = ', num2str(h2),' s)'])
ylabel('Displacement [m]')
legend('Newmark Integration Scheme')

subplot(4,1,4)
plot(t2,q(nodeA*Ndof-Ndof+2-18,:),'b',t,q2(nodeA*Ndof-Ndof+2-18,:),'g',t,q3(:,nodeA*Ndof-Ndof+2-18)','r')
title('Y Direction Comparisson')
xlabel(['Time [s] (Time increment = ', num2str(h2),' s)'])
ylabel('Displacement [m]')

%Displacement of the node A along z direction
figure
subplot(4,1,1)
plot(t,q(nodeA*Ndof-Ndof+3-18,:),'-b')
title('Dynamic Response - Z Direction')
xlabel(['Time [s] (Time increment = ', num2str(dt),' s)'])
ylabel('Displacement [m]')
legend('Modal Displacement Method')

subplot(4,1,2)
plot(t,q2(nodeA*Ndof-Ndof+3-18,:),'-g')
xlabel(['Time [s] (Time increment = ', num2str(dt),' s)'])
ylabel('Displacement [m]')
legend('Modal Acceleration Method')

subplot(4,1,3)
plot(t2,q3(:,nodeA*Ndof-Ndof+3-18)','-r')
xlabel(['Time [s] (Time increment = ', num2str(h2),' s)'])
ylabel('Displacement [m]')
legend('Newmark Integration Scheme')

subplot(4,1,4)
plot(t2,q(nodeA*Ndof-Ndof+3-18,:),'b',t,q2(nodeA*Ndof-Ndof+3-18,:),'g',t,q3(:,nodeA*Ndof-Ndof+3-18)','r')
title('Z Direction Comparisson')
xlabel(['Time [s] (Time increment = ', num2str(h2),' s)'])
ylabel('Displacement [m]')

%% Plotting of the responses in frequency domain
%%
xf = fft(q3(:,nodeA*Ndof-Ndof+1-18)'); %along x direction
xf(1) = [];
n = length(xf);
power = abs(xf(1:floor(n/2))).^2;
nyquist = 0.5;
fFFT = (1:n/2)/(n/2)*nyquist;
fFFT = fFFT*100;
figure()
subplot(3,1,1)
semilogy(fFFT,power,'b')
title('Fast Fouirer Transform - X direction')
xlabel('Frequency [Hz]')

yf = fft(q3(:,nodeA*Ndof-Ndof+2-18)'); %along y direction
yf(1) = [];
n = length(yf);
power = abs(yf(1:floor(n/2))).^2;
nyquist = 0.5;
fFFT = (1:n/2)/(n/2)*nyquist;
fFFT = fFFT*100;
subplot(3,1,2)
semilogy(fFFT,power,'b')
title('Fast Fouirer Transform - Y direction')
xlabel('Frequency [Hz]')

zf = fft(q3(:,nodeA*Ndof-Ndof+1-18)'); %along z direction
zf(1) = [];
n=length(zf);
power = abs(zf(1:floor(n/2))).^2;
nyquist = 0.5;
fFFT = (1:n/2)/(n/2)*nyquist;
fFFT = fFFT*100;
subplot(3,1,3)
semilogy(fFFT,power,'b')
title('Fast Fouirer Transform - Z direction')
xlabel('Frequency [Hz]')

%% MODEL REDUCTION
%%
Rnodes = [16,17,19,20]; %retained nodes
Rdofs  = [];            %retained degrees of freedom

for i=1:length(Rnodes)
    Rdofs=[Rdofs (Rnodes(i)*Ndof-Ndof+1-18):1:(Rnodes(i)*Ndof-18)];
end

Cdofs = 1:Ndofs;        %condensed degrees of freedom
Cdofs(Rdofs) = [];

kRR = K(Rdofs,Rdofs);   %retained part
kRC = K(Rdofs,Cdofs);          
kCC = K(Cdofs,Cdofs);   %condensed part
kCR = K(Cdofs,Rdofs);    

mRR = M(Rdofs,Rdofs);   %retained part
mRC = M(Rdofs,Cdofs);          
mCC = M(Cdofs,Cdofs);   %condensed part
mCR = M(Cdofs,Rdofs);    

%% Guyan-Iron Reduction
%%
Rgi = [eye(length(kRR)); -(kCC^-1)*kCR]; %transformation matrix

Kt=[kRR kRC;kCR kCC];
Mt=[mRR mRC;mCR mCC];
Kgi = Rgi'*Kt*Rgi;                       %reduced stiffness matrix
Mgi = Rgi'*Mt*Rgi;                       %reduced mass matrix

[Vgi, Dgi] = eigs(Kgi,Mgi,10,'smallestabs');
wGI = sqrt(diag(Dgi));                   %rad/s - first 10 angular eigenfrequencies
fGI = wGI/(2*pi);                        %Hz - first 10 eigenfrequencies

EigCompGI = zeros(1,length(fGI));        %eigenfrequency comparisson matrix for Guyan-Iron Reduction
for i=1:length(fGI)
    EigCompGI(i)=abs(fGI(i)-f(i))/f(i)*100;
end

%% Craig-Bampton Method
%%
[phi_r_CC, dCC] = eigs(kCC,mCC,10,'smallestabs');
wCC = sqrt(diag(dCC));
fCC = wCC/(2*pi); 

Neigenmodes = 4;
EigCompCB = zeros(length(f),Neigenmodes); %eigenfrequency comparisson matrix for Craig-Bampton Method
timeCB = zeros(Neigenmodes,1);
for i=1:Neigenmodes                       %accuracy vs. number of considered eigenmodes
    tic
    phi_r = phi_r_CC(:,1:i);

    Rcb2 = [zeros(length(kRR),i); phi_r]; %transformation submatrix
    Rcb = [Rgi Rcb2];                     %transformation matrix

    Kcb = Rcb'*Kt*Rcb;                    %reduced stiffness matrix
    Mcb = Rcb'*Mt*Rcb;                    %reduced mass matrix

    [Vcb, Dcb] = eigs(Kcb,Mcb,10,'smallestabs');
    wCB = sqrt(diag(Dcb));
    fCB = wCB/(2*pi);

    timeCB(i)=toc;
    
    for j=1:length(fCB)
        EigCompCB(j,i) = abs(fCB(j)-f(j))/f(j)*100;
    end

end

%% Plotting of the eigenfrequencies
%%
figure
plot(1:Neigenfrequency,f(1:Neigenfrequency),'-*','LineWidth', 2, 'MarkerSize', 10)
hold on
plot(1:Neigenfrequency,fGI(1:Neigenfrequency),'-d','LineWidth', 2, 'MarkerSize', 5)
hold on
plot(1:Neigenfrequency,fCB(1:Neigenfrequency),'-m*','LineWidth', 2, 'MarkerSize', 5)
title('Eigenfrequency comparison')
ylabel('frequency [Hz]')
xlabel('number of eigenfrequency')
legend('Full Model','Guyan-Iron','Craig Bampton (10 eigenmodes)')
grid on

figure
errorGI = (abs(fGI-f)./f)*100;
plot(1:1:Neigenfrequency,errorGI(1:Neigenfrequency),'-d','LineWidth', 2, 'MarkerSize', 5)
hold on
errorCB = (abs(fCB-f)./f)*100;
plot(1:1:Neigenfrequency,errorCB(1:Neigenfrequency),'-m*','LineWidth', 2, 'MarkerSize', 5)
title('Eigenfrequency error')
ylabel('error [%]')
xlabel('number of eigenfrequency')
legend('Guyan-Iron','Craig Bampton (4 eigenmodes)')
grid on
    
tEnd = toc(tStart)

%% FUNCTIONS
%%
function [I,J] = BeamProperties(A)

r = sqrt(A/pi);

I = pi*r^4/4;

J = I*2;

end

function L = BeamElementLength(x1,y1,z1,x2,y2,z2)

L = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);

end

function k = BeamElementStiffness(E,G,A,Iy,Iz,J,x1,y1,z1,x2,y2,z2,L)

w1 = E*A/L;
w2 = 12*E*Iz/(L*L*L);
w3 = 6*E*Iz/(L*L);
w4 = 4*E*Iz/L;
w5 = 2*E*Iz/L;
w6 = 12*E*Iy/(L*L*L);
w7 = 6*E*Iy/(L*L);
w8 = 4*E*Iy/L;
w9 = 2*E*Iy/L;
w10 = G*J/L;

kel = [w1 0 0 0 0 0 -w1 0 0 0 0 0 ;
        0 w2 0 0 0 w3 0 -w2 0 0 0 w3 ;
        0 0 w6 0 -w7 0 0 0 -w6 0 -w7 0 ;
        0 0 0 w10 0 0 0 0 0 -w10 0 0 ;
        0 0 -w7 0 w8 0 0 0 w7 0 w9 0 ;
        0 w3 0 0 0 w4 0 -w3 0 0 0 w5 ;
        -w1 0 0 0 0 0 w1 0 0 0 0 0 ;
        0 -w2 0 0 0 -w3 0 w2 0 0 0 -w3 ;
        0 0 -w6 0 w7 0 0 0 w6 0 w7 0 ;
        0 0 0 -w10 0 0 0 0 0 w10 0 0 ;
        0 0 -w7 0 w9 0 0 0 w7 0 w8 0 ;
        0 w3 0 0 0 w5 0 -w3 0 0 0 w4];

if x1 == x2 && y1 == y2
    if z2 > z1
        Lambda = [0 0 1 ; 0 1 0 ; -1 0 0];
    else
        Lambda = [0 0 -1 ; 0 1 0 ; 1 0 0];
    end
else
    CXx = (x2-x1)/L;
    CYx = (y2-y1)/L;
    CZx = (z2-z1)/L;
    D = sqrt(CXx*CXx + CYx*CYx);
    CXy = -CYx/D;
    CYy = CXx/D;
    CZy = 0;
    CXz = -CXx*CZx/D;
    CYz = -CYx*CZx/D;
    CZz = D;
    Lambda = [CXx CYx CZx ; CXy CYy CZy ; CXz CYz CZz];
end

R = [Lambda zeros(3) zeros(3) zeros(3) ;
    zeros(3) Lambda zeros(3) zeros(3) ;
    zeros(3) zeros(3) Lambda zeros(3) ;
    zeros(3) zeros(3) zeros(3) Lambda];

k = R'*kel*R;

end

function m = BeamElementMassMatrix(rho,A,x1,y1,z1,x2,y2,z2,L,r)
   
mel = rho*A*L*...
        [1/3 0 0 0 0 0 1/6 0 0 0 0 0;
        0 13/35 0 0 0 11*L/210 0 9/70 0 0 0 -13*L/420;
        0 0 13/35 0 -11*L/210 0 0 0 9/70 0 13*L/420 0;
        0 0 0 r^2/3 0 0 0 0 0 r^2/6 0 0;
        0 0 -11*L/210 0 L^2/105 0 0 0 -13*L/420 0 -L^2/140 0;
        0 11*L/210 0 0 0 L^2/105 0 13*L/420 0 0 0 -L^2/140;
        1/6 0 0 0 0 0 1/3 0 0 0 0 0;
        0 9/70 0 0 0 13*L/420 0 13/35 0 0 0 -11*L/210;
        0 0 9/70 0 -13*L/420 0 0 0 13/35 0 11*L/210 0;
        0 0 0 r^2/6 0 0 0 0 0 r^2/3 0 0;
        0 0 13*L/420 0 -L^2/140 0 0 0 11*L/210 0 L^2/105 0;
        0 -13*L/420 0 0 0 -L^2/140 0 -11*L/210 0 0 0 L^2/105];

if x1 == x2 && y1 == y2
    if z2 > z1
        Lambda = [0 0 1 ; 0 1 0 ; -1 0 0];
    else
        Lambda = [0 0 -1 ; 0 1 0 ; 1 0 0];
    end
else
    CXx = (x2-x1)/L;
    CYx = (y2-y1)/L;
    CZx = (z2-z1)/L;
    D = sqrt(CXx*CXx + CYx*CYx);
    CXy = -CYx/D;
    CYy = CXx/D;
    CZy = 0;
    CXz = -CXx*CZx/D;
    CYz = -CYx*CZx/D;
    CZz = D;
    Lambda = [CXx CYx CZx ; CXy CYy CZy ; CXz CYz CZz];
end

R = [Lambda zeros(3) zeros(3) zeros(3) ;
    zeros(3) Lambda zeros(3) zeros(3) ;
    zeros(3) zeros(3) Lambda zeros(3) ;
    zeros(3) zeros(3) zeros(3) Lambda];

m = R'*mel*R;

end

function K = CraneJibAssemble(k,K,i,j)

K(6*i-5,6*i-5) = K(6*i-5,6*i-5) + k(1,1);
K(6*i-5,6*i-4) = K(6*i-5,6*i-4) + k(1,2);
K(6*i-5,6*i-3) = K(6*i-5,6*i-3) + k(1,3);
K(6*i-5,6*i-2) = K(6*i-5,6*i-2) + k(1,4);
K(6*i-5,6*i-1) = K(6*i-5,6*i-1) + k(1,5);
K(6*i-5,6*i) = K(6*i-5,6*i) + k(1,6);
K(6*i-5,6*j-5) = K(6*i-5,6*j-5) + k(1,7);
K(6*i-5,6*j-4) = K(6*i-5,6*j-4) + k(1,8);
K(6*i-5,6*j-3) = K(6*i-5,6*j-3) + k(1,9);
K(6*i-5,6*j-2) = K(6*i-5,6*j-2) + k(1,10);
K(6*i-5,6*j-1) = K(6*i-5,6*j-1) + k(1,11);
K(6*i-5,6*j) = K(6*i-5,6*j) + k(1,12);
K(6*i-4,6*i-5) = K(6*i-4,6*i-5) + k(2,1);
K(6*i-4,6*i-4) = K(6*i-4,6*i-4) + k(2,2);
K(6*i-4,6*i-3) = K(6*i-4,6*i-3) + k(2,3);
K(6*i-4,6*i-2) = K(6*i-4,6*i-2) + k(2,4);
K(6*i-4,6*i-1) = K(6*i-4,6*i-1) + k(2,5);
K(6*i-4,6*i) = K(6*i-4,6*i) + k(2,6);
K(6*i-4,6*j-5) = K(6*i-4,6*j-5) + k(2,7);
K(6*i-4,6*j-4) = K(6*i-4,6*j-4) + k(2,8);
K(6*i-4,6*j-3) = K(6*i-4,6*j-3) + k(2,9);
K(6*i-4,6*j-2) = K(6*i-4,6*j-2) + k(2,10);
K(6*i-4,6*j-1) = K(6*i-4,6*j-1) + k(2,11);
K(6*i-4,6*j) = K(6*i-4,6*j) + k(2,12);
K(6*i-3,6*i-5) = K(6*i-3,6*i-5) + k(3,1);
K(6*i-3,6*i-4) = K(6*i-3,6*i-4) + k(3,2);
K(6*i-3,6*i-3) = K(6*i-3,6*i-3) + k(3,3);
K(6*i-3,6*i-2) = K(6*i-3,6*i-2) + k(3,4);
K(6*i-3,6*i-1) = K(6*i-3,6*i-1) + k(3,5);
K(6*i-3,6*i) = K(6*i-3,6*i) + k(3,6);
K(6*i-3,6*j-5) = K(6*i-3,6*j-5) + k(3,7);
K(6*i-3,6*j-4) = K(6*i-3,6*j-4) + k(3,8);
K(6*i-3,6*j-3) = K(6*i-3,6*j-3) + k(3,9);
K(6*i-3,6*j-2) = K(6*i-3,6*j-2) + k(3,10);
K(6*i-3,6*j-1) = K(6*i-3,6*j-1) + k(3,11);
K(6*i-3,6*j) = K(6*i-3,6*j) + k(3,12);
K(6*i-2,6*i-5) = K(6*i-2,6*i-5) + k(4,1);
K(6*i-2,6*i-4) = K(6*i-2,6*i-4) + k(4,2);
K(6*i-2,6*i-3) = K(6*i-2,6*i-3) + k(4,3);
K(6*i-2,6*i-2) = K(6*i-2,6*i-2) + k(4,4);
K(6*i-2,6*i-1) = K(6*i-2,6*i-1) + k(4,5);
K(6*i-2,6*i) = K(6*i-2,6*i) + k(4,6);
K(6*i-2,6*j-5) = K(6*i-2,6*j-5) + k(4,7);
K(6*i-2,6*j-4) = K(6*i-2,6*j-4) + k(4,8);
K(6*i-2,6*j-3) = K(6*i-2,6*j-3) + k(4,9);
K(6*i-2,6*j-2) = K(6*i-2,6*j-2) + k(4,10);
K(6*i-2,6*j-1) = K(6*i-2,6*j-1) + k(4,11);
K(6*i-2,6*j) = K(6*i-2,6*j) + k(4,12);
K(6*i-1,6*i-5) = K(6*i-1,6*i-5) + k(5,1);
K(6*i-1,6*i-4) = K(6*i-1,6*i-4) + k(5,2);
K(6*i-1,6*i-3) = K(6*i-1,6*i-3) + k(5,3);
K(6*i-1,6*i-2) = K(6*i-1,6*i-2) + k(5,4);
K(6*i-1,6*i-1) = K(6*i-1,6*i-1) + k(5,5);
K(6*i-1,6*i) = K(6*i-1,6*i) + k(5,6);
K(6*i-1,6*j-5) = K(6*i-1,6*j-5) + k(5,7);
K(6*i-1,6*j-4) = K(6*i-1,6*j-4) + k(5,8);
K(6*i-1,6*j-3) = K(6*i-1,6*j-3) + k(5,9);
K(6*i-1,6*j-2) = K(6*i-1,6*j-2) + k(5,10);
K(6*i-1,6*j-1) = K(6*i-1,6*j-1) + k(5,11);
K(6*i-1,6*j) = K(6*i-1,6*j) + k(5,12);
K(6*i,6*i-5) = K(6*i,6*i-5) + k(6,1);
K(6*i,6*i-4) = K(6*i,6*i-4) + k(6,2);
K(6*i,6*i-3) = K(6*i,6*i-3) + k(6,3);
K(6*i,6*i-2) = K(6*i,6*i-2) + k(6,4);
K(6*i,6*i-1) = K(6*i,6*i-1) + k(6,5);
K(6*i,6*i) = K(6*i,6*i) + k(6,6);
K(6*i,6*j-5) = K(6*i,6*j-5) + k(6,7);
K(6*i,6*j-4) = K(6*i,6*j-4) + k(6,8);
K(6*i,6*j-3) = K(6*i,6*j-3) + k(6,9);
K(6*i,6*j-2) = K(6*i,6*j-2) + k(6,10);
K(6*i,6*j-1) = K(6*i,6*j-1) + k(6,11);
K(6*i,6*j) = K(6*i,6*j) + k(6,12);
K(6*j-5,6*i-5) = K(6*j-5,6*i-5) + k(7,1);
K(6*j-5,6*i-4) = K(6*j-5,6*i-4) + k(7,2);
K(6*j-5,6*i-3) = K(6*j-5,6*i-3) + k(7,3);
K(6*j-5,6*i-2) = K(6*j-5,6*i-2) + k(7,4);
K(6*j-5,6*i-1) = K(6*j-5,6*i-1) + k(7,5);
K(6*j-5,6*i) = K(6*j-5,6*i) + k(7,6);
K(6*j-5,6*j-5) = K(6*j-5,6*j-5) + k(7,7);
K(6*j-5,6*j-4) = K(6*j-5,6*j-4) + k(7,8);
K(6*j-5,6*j-3) = K(6*j-5,6*j-3) + k(7,9);
K(6*j-5,6*j-2) = K(6*j-5,6*j-2) + k(7,10);
K(6*j-5,6*j-1) = K(6*j-5,6*j-1) + k(7,11);
K(6*j-5,6*j) = K(6*j-5,6*j) + k(7,12);
K(6*j-4,6*i-5) = K(6*j-4,6*i-5) + k(8,1);
K(6*j-4,6*i-4) = K(6*j-4,6*i-4) + k(8,2);
K(6*j-4,6*i-3) = K(6*j-4,6*i-3) + k(8,3);
K(6*j-4,6*i-2) = K(6*j-4,6*i-2) + k(8,4);
K(6*j-4,6*i-1) = K(6*j-4,6*i-1) + k(8,5);
K(6*j-4,6*i) = K(6*j-4,6*i) + k(8,6);
K(6*j-4,6*j-5) = K(6*j-4,6*j-5) + k(8,7);
K(6*j-4,6*j-4) = K(6*j-4,6*j-4) + k(8,8);
K(6*j-4,6*j-3) = K(6*j-4,6*j-3) + k(8,9);
K(6*j-4,6*j-2) = K(6*j-4,6*j-2) + k(8,10);
K(6*j-4,6*j-1) = K(6*j-4,6*j-1) + k(8,11);
K(6*j-4,6*j) = K(6*j-4,6*j) + k(8,12);
K(6*j-3,6*i-5) = K(6*j-3,6*i-5) + k(9,1);
K(6*j-3,6*i-4) = K(6*j-3,6*i-4) + k(9,2);
K(6*j-3,6*i-3) = K(6*j-3,6*i-3) + k(9,3);
K(6*j-3,6*i-2) = K(6*j-3,6*i-2) + k(9,4);
K(6*j-3,6*i-1) = K(6*j-3,6*i-1) + k(9,5);
K(6*j-3,6*i) = K(6*j-3,6*i) + k(9,6);
K(6*j-3,6*j-5) = K(6*j-3,6*j-5) + k(9,7);
K(6*j-3,6*j-4) = K(6*j-3,6*j-4) + k(9,8);
K(6*j-3,6*j-3) = K(6*j-3,6*j-3) + k(9,9);
K(6*j-3,6*j-2) = K(6*j-3,6*j-2) + k(9,10);
K(6*j-3,6*j-1) = K(6*j-3,6*j-1) + k(9,11);
K(6*j-3,6*j) = K(6*j-3,6*j) + k(9,12);
K(6*j-2,6*i-5) = K(6*j-2,6*i-5) + k(10,1);
K(6*j-2,6*i-4) = K(6*j-2,6*i-4) + k(10,2);
K(6*j-2,6*i-3) = K(6*j-2,6*i-3) + k(10,3);
K(6*j-2,6*i-2) = K(6*j-2,6*i-2) + k(10,4);
K(6*j-2,6*i-1) = K(6*j-2,6*i-1) + k(10,5);
K(6*j-2,6*i) = K(6*j-2,6*i) + k(10,6);
K(6*j-2,6*j-5) = K(6*j-2,6*j-5) + k(10,7);
K(6*j-2,6*j-4) = K(6*j-2,6*j-4) + k(10,8);
K(6*j-2,6*j-3) = K(6*j-2,6*j-3) + k(10,9);
K(6*j-2,6*j-2) = K(6*j-2,6*j-2) + k(10,10);
K(6*j-2,6*j-1) = K(6*j-2,6*j-1) + k(10,11);
K(6*j-2,6*j) = K(6*j-2,6*j) + k(10,12);
K(6*j-1,6*i-5) = K(6*j-1,6*i-5) + k(11,1);
K(6*j-1,6*i-4) = K(6*j-1,6*i-4) + k(11,2);
K(6*j-1,6*i-3) = K(6*j-1,6*i-3) + k(11,3);
K(6*j-1,6*i-2) = K(6*j-1,6*i-2) + k(11,4);
K(6*j-1,6*i-1) = K(6*j-1,6*i-1) + k(11,5);
K(6*j-1,6*i) = K(6*j-1,6*i) + k(11,6);
K(6*j-1,6*j-5) = K(6*j-1,6*j-5) + k(11,7);
K(6*j-1,6*j-4) = K(6*j-1,6*j-4) + k(11,8);
K(6*j-1,6*j-3) = K(6*j-1,6*j-3) + k(11,9);
K(6*j-1,6*j-2) = K(6*j-1,6*j-2) + k(11,10);
K(6*j-1,6*j-1) = K(6*j-1,6*j-1) + k(11,11);
K(6*j-1,6*j) = K(6*j-1,6*j) + k(11,12);
K(6*j,6*i-5) = K(6*j,6*i-5) + k(12,1);
K(6*j,6*i-4) = K(6*j,6*i-4) + k(12,2);
K(6*j,6*i-3) = K(6*j,6*i-3) + k(12,3);
K(6*j,6*i-2) = K(6*j,6*i-2) + k(12,4);
K(6*j,6*i-1) = K(6*j,6*i-1) + k(12,5);
K(6*j,6*i) = K(6*j,6*i) + k(12,6);
K(6*j,6*j-5) = K(6*j,6*j-5) + k(12,7);
K(6*j,6*j-4) = K(6*j,6*j-4) + k(12,8);
K(6*j,6*j-3) = K(6*j,6*j-3) + k(12,9);
K(6*j,6*j-2) = K(6*j,6*j-2) + k(12,10);
K(6*j,6*j-1) = K(6*j,6*j-1) + k(12,11);
K(6*j,6*j) = K(6*j,6*j) + k(12,12);

end
