%Generalized Patch Dynamics Scheme of type-I in Equation-free Multiscale Modelling%
     % Author: Tanay Kumar Karmakar & Durga Charan Dalal
     % Indian Institute of Technology Guwahati, India.
     % Date: April 9, 2024


%% Analytical Solution

%% Domain and time 
La=1;       % Length of the domain   
Ta=1;       % Final time 
Nta=100;    % Number of macro time levels       
Nxa=20;     % Number of macro grids 

%% The IC: 
xa=linspace(0,La,Nxa+1); % Space discretization
ta=linspace(0,Ta,Nta+1); % Time discretization
ua=zeros(Nxa+1,Nta+1);   % Initial declaration

for ka=1:Nta+1           
    for ia=1:Nxa+1       
        ua(ia,ka)=(sin(pi*xa(ia)))*exp((1-pi*pi)*ta(ka));  
    end
end

figure(1)
plot(xa,ua(:,round(Nta/100)+1),'--^c',xa,ua(:,round(Nta/10)+1),'--^r',xa,ua(:,round(Nta/5)+1),'--^g');
xlabel('x')
ylabel('Solution')
grid on
hold on;

figure(2)
plot(ta,ua(1,:),'--c',ta,ua(round(Nxa/4)+1,:),'--r',ta,ua(round(Nxa/2)+1,:),'--g');
xlabel('t')
ylabel('Solution')
grid on
hold on;



tic
%% Solving 1D linear reaction-diffusion equation by PDS

%% Macro domain and time, it's discretization
Li=La;          
Ti=Ta;         

Nt=Nta;       
Dt=Ti/Nt;       
Nx=Nxa;         
Dx=Li/Nx;      
R=Dt/(Dx*Dx);     
X=linspace(0,Li,Nx+1);
T=linspace(0,Ti,Nt+1);
disp(R);

%% Patch size and it's discretization
h=0.02;         % Patch size
tau=1e-5;       % Time of micro simulation
nt=500; %20000; % #micro time levels
dt=tau/nt;      % Micro time level size
nx=20; %600;    % #micro grids
dx=h/nx;        % Length of micro grids
r=dt/(dx*dx);   % Stability criterion
disp(r);

%% Distribution of GTTs
% Uniform GTTs
%ktau=[14,14];
%ktau=[10,10,10];
%ktau=[8,8,8,8];
%ktau=[7,7,7,7,6];
%ktau=[6,6,6,6,6,6];
ktau=[5,5,5,5,5,5,5,5];
%ktau=[4,4,4,4,4,4,4,4,4,4,4,4];             
%ktau=[3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3];

% Non-uniform GTTs
%ktau=[6,4,10,3,7,6];

%% Distribution of extrapolation time
% Uniform extrapolations
E_tot=Dt-24*tau;
Ng=numel(ktau); % #GTS groups
Ex=E_tot/Ng;


% % Non-uniform extrapolations
% E_tot=Dt-24*tau;
% Ex=[E_tot*0.1, E_tot*0.15, E_tot*0.05, E_tot*0.25, E_tot*0.05, E_tot*0.4];
% Ng=numel(ktau); % #GTS groups
% Ne=numel(Ex);


%% The IC at macro level
for j=1:Nx+1
    U(j,1,1,1)=sin(pi*X(j)); 
end

%% The BCs at macro level
for i=1:Nt+1              
    for l=1:Ng
        for k=1:ktau(l)
            U(1,i,k,l)=0;
            U(Nx+1,i,k,l)=0;  
        end
    end 
end




%% Loop
for i=2:Nt+1        
  for l=1:Ng
    for k=2:ktau(l)
        for j=2:Nx  

        
        % PtBCs
             sl(j)=((Dx-h)*U(j+1,i-1,k-1,l)+2*h*U(j,i-1,k-1,l)-(Dx+h)*U(j-1,i-1,k-1,l))/(2*Dx*Dx); % Left BC
             sr(j)=((Dx+h)*U(j+1,i-1,k-1,l)-2*h*U(j,i-1,k-1,l)-(Dx-h)*U(j-1,i-1,k-1,l))/(2*Dx*Dx); % Right BC
             

        % Lifting
             for q=1:nx+1    
                 x(q)=X(j)-h/2+(q-1)*dx; 
                 u(q,1)=((U(j+1,i-1,k-1,l)-2*U(j,i-1,k-1,l)+U(j-1,i-1,k-1,l))*(x(q)-X(j))^2)/(2*Dx*Dx)+((U(j+1,i-1,k-1,l)-U(j-1,i-1,k-1,l))*(x(q)-X(j)))/(2*Dx)+U(j,i-1,k-1,l)-h*h*(U(j+1,i-1,k-1,l)-2*U(j,i-1,k-1,l)+U(j-1,i-1,k-1,l))/(24*Dx*Dx);
                 
             end
        
           
        % FTCS (Evolving)
             for p=2:nt+1       % Time loop
                 for q=2:nx       % Space loop
                     u(q,p)=u(q,p-1)+r*(u(q-1,p-1)+u(q+1,p-1)-2*u(q,p-1))+dt*u(q,p-1);  
                 end    
                 u(1,p)=u(2,p)-sl(j)*dx;
                 u(nx+1,p)=u(nx,p)+sr(j)*dx;
             end
             
             
             
        % Restricting
             % Simpson 1/3 
             v=0;
             for q=2:2:nx
                 v=v+4*u(q,nt+1);
             end
             for q=3:2:nx-1
                 v=v+2*u(q,nt+1);
             end
             v=v+u(1,nt+1)+u(nx+1,nt+1);
             U(j,i-1,k,l)=dx*v/(3*h);
             


        end    

    end        

    %% Extrapolation
     % Uniform extrapolations
     for j=2:Nx
         U(j,i-1,1,l+1)=U(j,i-1,ktau(l)-1,l)+Ex*(U(j,i-1,ktau(l),l)-U(j,i-1,ktau(l)-1,l))/(tau); 
     end

     % % Non-uniform extrapolations
     % for j=2:Nx
     %     U(j,i-1,1,l+1)=U(j,i-1,ktau(l)-1,l)+Ex(l)*(U(j,i-1,ktau(l),l)-U(j,i-1,ktau(l)-1,l))/(tau); 
     % end

  end     

      

    % Final coarse level value
     for j=2:Nx
         U(j,i,1,1)=U(j,i-1,1,Ng+1);
     end

    

   
end         

toc

M_Per=max(max(abs(U(2:Nx,:,1,1)-ua(2:Nx,:))./ua(2:Nx,:)))*100;
fprintf('Maximum percentage error=%4.6f\n', M_Per);



    figure(1)
    plot(X,U(:,1,1,1),'-*r', X,U(:,round(Nt/100)+1,1,1),'-*g',X,U(:,round(Nt/10)+1,1,1),'-*b',X,U(:,round(Nt/5)+1,1,1),'-*m',X,U(:,Nt+1,1,1),'-*k');
    hold on;

    figure(2)
    plot(T,U(1,:,1,1),'-b',T,U(round(Nx/4)+1,:,1,1),'-m',T,U(round(Nx/2)+1,:,1,1),'-k');
    hold on;

   figure(3)
    contourf(U(:,:,1,1),'ShowText','on');

    figure(4)
    contourf(abs(U(:,:,1,1)-ua(:,:)),'ShowText','on');

    figure(5)
    contourf((abs(U(2:Nx,:,1,1)-ua(2:Nx,:))./ua(2:Nx,:))*100,'ShowText','on');

figure(6)
tan1=U(:,:,1,1);
grid on
colormap jet
colorbar
mesh(X,T,tan1')
title('Full Error New');
xlabel('x')
ylabel('t')
zlabel('E')

figure(7)
tan2=abs(U(:,:,1,1)-ua(:,:));
grid on
colormap jet
colorbar
mesh(X,T,tan2')
title('Full Error New');
xlabel('x')
ylabel('t')
zlabel('E')

figure(8)
tan3=(abs(U(2:Nx,:,1,1)-ua(2:Nx,:))./ua(2:Nx,:));
grid on
colormap jet
colorbar
mesh(X(2:Nx),T,tan3')
title('Full Error New');
xlabel('x')
ylabel('t')
zlabel('E')


figure(9)
surf(X,T,U(:,:,1,1)','FaceColor','interp',...
'EdgeColor','none',...
'FaceLighting','gouraud');
%daspect([5 5 1])
axis tight
view(-50,30)
camlight right

figure(10)
surf(X,T,abs(U(:,:,1,1)-ua(:,:))','FaceColor','interp',...
'EdgeColor','none',...
'FaceLighting','gouraud');
%daspect([5 5 1])
axis tight
view(-50,30)
camlight right

figure(11)
surf(X(2:Nx),T,(abs(U(2:Nx,:,1,1)-ua(2:Nx,:))./ua(2:Nx,:))'*100,'FaceColor','interp',...
'EdgeColor','none',...
'FaceLighting','gouraud');
%daspect([5 5 1])
axis tight
view(-50,30)
camlight right
          
