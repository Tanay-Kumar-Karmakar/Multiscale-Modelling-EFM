   %--Usual Patch Dynamics Scheme in Equation-free Multiscale Modelling--%
     % Author: Tanay Kumar Karmakar & Durga Charan Dalal
     % Indian Institute of Technology Guwahati, India.
     % Date: April 9, 2024


%% Analytical Solution

%% Domain and time 
La=1;          % Length of the domain
Ta=1;          % Final time
Nta=1000;      % Number of macro time levels  
Nxa=20;        % Number of macro grids

%% The IC:
xa=linspace(0,La,Nxa+1);  % Space discretization
ta=linspace(0,Ta,Nta+1);  % Time discretization
ua=zeros(Nxa+1,Nta+1);    % Initial declaration

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
h=0.02;        % Patch size
tau=1e-5;      % Time of micro computation
ktau=26;       % Number of GTTs
nt=20000;      % #micro time levels
dt=tau/nt;     % Micro time level size
nx=600;        % #micro grids
dx=h/nx;       % Length of micro grids
r=dt/(dx*dx);  % Stability criterion
disp(r);

U=zeros(Nx+1,Nt+1,ktau);  % Initial declaration
%% The IC at macro level
for j=1:Nx+1
    U(j,1,1)=sin(pi*X(j)); 
end

%% The BCs at macro level
for i=1:Nt+1
    for k=1:ktau
        U(1,i,k)=0;
        U(Nx+1,i,k)=0;
    end
end


%% Loop
for i=2:Nt+1        
    for k=2:ktau
        for j=2:Nx  
        
        % PtBCs
             sl(j)=((Dx-h)*U(j+1,i-1,k-1)+2*h*U(j,i-1,k-1)-(Dx+h)*U(j-1,i-1,k-1))/(2*Dx*Dx); % Left BC for the patch
             sr(j)=((Dx+h)*U(j+1,i-1,k-1)-2*h*U(j,i-1,k-1)-(Dx-h)*U(j-1,i-1,k-1))/(2*Dx*Dx); % Right BC for the patch
            

        % Lifting
             for q=1:nx+1    
                 x(q)=X(j)-h/2+(q-1)*dx; 
                 u(q,1)=((U(j+1,i-1,k-1)-2*U(j,i-1,k-1)+U(j-1,i-1,k-1))*(x(q)-X(j))^2)/(2*Dx*Dx)+((U(j+1,i-1,k-1)-U(j-1,i-1,k-1))*(x(q)-X(j)))/(2*Dx)+U(j,i-1,k-1)-h*h*(U(j+1,i-1,k-1)-2*U(j,i-1,k-1)+U(j-1,i-1,k-1))/(24*Dx*Dx);
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
             
             % Simpson's 1/3 
             v=0;
             for q=2:2:nx
                 v=v+4*u(q,nt+1);
             end
             for q=3:2:nx
                 v=v+2*u(q,nt+1);
             end
             v=v+u(1,nt+1)+u(nx+1,nt+1);
             U(j,i-1,k)=dx*v/(3*h);
         
        end   

    end        

    %% Extrapolation
     for j=2:Nx
         U(j,i,1)=U(j,i-1,ktau-1)+(Dt-(ktau-2)*tau)*(U(j,i-1,ktau)-U(j,i-1,ktau-1))/(tau); 
     end
   
end          
   

 M_Per=max(max(abs(U(2:Nx,:,1)-ua(2:Nx,:))./ua(2:Nx,:)))*100;
fprintf('Maximum percentage error=%4.6f\n', M_Per);


figure(1)
plot(X,U(:,1,1),'-pr',X,U(:,round(Nt/100)+1,1),'-pg',X,U(:,round(Nt/10)+1,1),'-pb',X,U(:,round(Nt/5)+1,1),'-pm',X,U(:,Nt+1,1),'-pk');
hold on;

figure(2)
plot(T,U(1,:,1),'-g',T,U(round(Nx/4)+1,:,1),'-b',T,U(round(Nx/2)+1,:,1),'-m');
hold on;



toc