% Program that implements the 5th order method for stiff problems called
% expRK5s10 that is constructed in 
% V. T. Luan, Efficient exponential Runge-Kutta methods of high order:
% construction and implementation, BIT Numerical Mathematics,61 (2021)
% 535-560
% 
% without avoiding the order reduction 
% u_t=u_xx+u_yy+g(t,u), (x,y)\in [0,1]x[0,1] with boundary conditions
% u(t,0,y)=gx0(t,y), u(t,1,y)=gx1(t,y), u(t,x,0)=gy0(t,x), and u(t,x,1)=gy1(t,x) 
% with g(t,u)=u^2+h(t,u) and with exact solution u(t,x,y)=cos(t+x+y)
% For the spatial discretization we have used the 9 point formula. 
%

% Some initial values, such as h and some more that are used several times 
% along the program
% JJ is such that h=1/JJ is the grid diameter in [0,1] for x and y

JJ=80;
dim=JJ-1;
dim2=(JJ-1)*(JJ-1);
h=1/JJ;
Chtilde=zeros(dim2,1);
Dhtilde=zeros(dim2,1);
K2=zeros(dim2,1);
K3=zeros(dim2,1);
K4=zeros(dim2,1);
K5=zeros(dim2,1);
K6=zeros(dim2,1);
K7=zeros(dim2,1);
K8=zeros(dim2,1);
vech=zeros(dim2,1);
vechk2=zeros(dim2,1);
vechk=zeros(dim2,1);
vechk4=zeros(dim2,1);
vechk5=zeros(dim2,1);
vechk310=zeros(dim2,1); 
ERROR1=zeros(dim,dim);
funu=zeros(dim2,1);
vecb2=zeros(dim2,2);
vecb3=zeros(dim2,3);
vecb4=zeros(dim2,4);   
vecb5=zeros(dim2,5);

x=zeros(dim,1);
y=zeros(dim,1);
x=[h:h:1]';
y=x;

% Matrix A_{h,0} is given by M^{-1}A. M^{-1}, is not computed, when
% necessary, a system is solved, phipmM9points and phipm_simul_iom9points 
% are used. Matrices A and M are not computed because they are very large.

% n is such that the time step size is k=1/n.
n=2;
k=1/n;
Tf=n*k;
h2=h^2;
mult=12/h2;
mult2=h2/12;

% The program runs for 8 different values of k, from k=1/2, in order to calculate
% the error and the order of the method
for ll=1:8
    
    k=1/n;
    k2=k/2;
    k3=k/3;
    k4=k/4;
    k310=3*k/10;
    k34=3*k/4;
    
    t=0;       
       
    % U is the initial value of u(0,x,y)
    U=zeros(dim2,1);
    for ii=1:dim
        for jj=1:dim
            U(ii+(jj-1)*dim,1)=cos(x(ii)+y(jj));
        end
    end
    
    % CPU time calculus starts
    tstart=tic;
   
    % For the local error r=1. For the global one, r=n         
    for kk=1:n
           
        bb=2/3;
        cc=1/6;  
        
        % Ch and Dh are calculated by solving a system M x=u. Here, u is
        % calculated u and then, system M x=u is solved. 
        % There are two types of boundary values, Ch g and D_h g. 
        % Boundaries of the form C_h g are of the form (12/h^2) M^{-1} bound Ag 
        % and the ones of the form D_h g are M^{-1} fron M g. 
        
        % All the vector are (JJ-1)*(JJ-1). The result corresponding to
        % (x(i),y(JJ)) is at position ii+JJ(jj-1)
                      
        % Boundary Chu
        MM1=cos(t+y(1))+cos(t+x(1));
        MM2=cos(t+x(2))+cos(t+y(2))+cos(t);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+x(ii));
            MM2=cos(t+x(ii+1))+cos(t+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+1+y(1))+cos(t+x(dim));
        MM2=cos(t+1+y(2))+cos(t+1)+cos(t+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;
           
        for m=2:dim-1
            MM1=cos(t+y(m));
            MM2=cos(t+y(m+1))+cos(t+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            
            MM1=cos(t+1+y(m));
            MM2=cos(t+1+y(m+1))+cos(t+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end
           
        MM1=cos(t+x(1)+1)+cos(t+y(dim));
        MM2=cos(t+x(2)+1)+cos(t+1)+cos(t+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+x(ii)+1);
            MM2=cos(t+x(ii+1)+1)+cos(t+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+1+y(dim))+cos(t+x(dim)+1);
        MM2=cos(t+2)+cos(t+1+y(dim-1))+cos(t+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        
        Chu=mult*multMm1(dim,Chtilde); 
                
        % Boundary Chu(tn+k2)
        MM1=cos(t+k2+y(1))+cos(t+k2+x(1));
        MM2=cos(t+k2+x(2))+cos(t+k2+y(2))+cos(t+k2);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k2+x(ii));
            MM2=cos(t+k2+x(ii+1))+cos(t+k2+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k2+1+y(1))+cos(t+k2+x(dim));
        MM2=cos(t+k2+1+y(2))+cos(t+k2+1)+cos(t+k2+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;
           
        for m=2:dim-1
            MM1=cos(t+k2+y(m));
            MM2=cos(t+k2+y(m+1))+cos(t+k2+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            
            MM1=cos(t+k2+1+y(m));
            MM2=cos(t+k2+1+y(m+1))+cos(t+k2+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end
           
        MM1=cos(t+k2+x(1)+1)+cos(t+k2+y(dim));
        MM2=cos(t+k2+x(2)+1)+cos(t+k2+1)+cos(t+k2+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k2+x(ii)+1);
            MM2=cos(t+k2+x(ii+1)+1)+cos(t+k2+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k2+1+y(dim))+cos(t+k2+x(dim)+1);
        MM2=cos(t+k2+2)+cos(t+k2+1+y(dim-1))+cos(t+k2+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        
        Chuk2=mult*multMm1(dim,Chtilde); 
                
        % Boundary Chu(tn+k)
        MM1=cos(t+k+y(1))+cos(t+k+x(1));
        MM2=cos(t+k+x(2))+cos(t+k+y(2))+cos(t+k);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k+x(ii));
            MM2=cos(t+k+x(ii+1))+cos(t+k+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k+1+y(1))+cos(t+k+x(dim));
        MM2=cos(t+k+1+y(2))+cos(t+k+1)+cos(t+k+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;           
        for m=2:dim-1
            MM1=cos(t+k+y(m));
            MM2=cos(t+k+y(m+1))+cos(t+k+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            
            MM1=cos(t+k+1+y(m));
            MM2=cos(t+k+1+y(m+1))+cos(t+k+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end           
        MM1=cos(t+k+x(1)+1)+cos(t+k+y(dim));
        MM2=cos(t+k+x(2)+1)+cos(t+k+1)+cos(t+k+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k+x(ii)+1);
            MM2=cos(t+k+x(ii+1)+1)+cos(t+k+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k+1+y(dim))+cos(t+k+x(dim)+1);
        MM2=cos(t+k+2)+cos(t+k+1+y(dim-1))+cos(t+k+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        
        Chuk=mult*multMm1(dim,Chtilde); 
                           
        % Boundary Chu(tn+k3)
        MM1=cos(t+k3+y(1))+cos(t+k3+x(1));
        MM2=cos(t+k3+x(2))+cos(t+k3+y(2))+cos(t+k3);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k3+x(ii));
            MM2=cos(t+k3+x(ii+1))+cos(t+k3+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k3+1+y(1))+cos(t+k3+x(dim));
        MM2=cos(t+k3+1+y(2))+cos(t+k3+1)+cos(t+k3+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;           
        for m=2:dim-1
            MM1=cos(t+k3+y(m));
            MM2=cos(t+k3+y(m+1))+cos(t+k3+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            
            MM1=cos(t+k3+1+y(m));
            MM2=cos(t+k3+1+y(m+1))+cos(t+k3+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end           
        MM1=cos(t+k3+x(1)+1)+cos(t+k3+y(dim));
        MM2=cos(t+k3+x(2)+1)+cos(t+k3+1)+cos(t+k3+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k3+x(ii)+1);
            MM2=cos(t+k3+x(ii+1)+1)+cos(t+k3+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k3+1+y(dim))+cos(t+k3+x(dim)+1);
        MM2=cos(t+k3+2)+cos(t+k3+1+y(dim-1))+cos(t+k3+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        
        Chuk3=mult*multMm1(dim,Chtilde); 
        
        
        % Boundary Chu(tn+k4)        
        MM1=cos(t+k4+y(1))+cos(t+k4+x(1));
        MM2=cos(t+k4+x(2))+cos(t+k4+y(2))+cos(t+k4);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k4+x(ii));
            MM2=cos(t+k4+x(ii+1))+cos(t+k4+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k4+1+y(1))+cos(t+k4+x(dim));
        MM2=cos(t+k4+1+y(2))+cos(t+k4+1)+cos(t+k4+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;           
        for m=2:dim-1
            MM1=cos(t+k4+y(m));
            MM2=cos(t+k4+y(m+1))+cos(t+k4+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;            
            MM1=cos(t+k4+1+y(m));
            MM2=cos(t+k4+1+y(m+1))+cos(t+k4+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end           
        MM1=cos(t+k4+x(1)+1)+cos(t+k4+y(dim));
        MM2=cos(t+k4+x(2)+1)+cos(t+k4+1)+cos(t+k4+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k4+x(ii)+1);
            MM2=cos(t+k4+x(ii+1)+1)+cos(t+k4+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k4+1+y(dim))+cos(t+k4+x(dim)+1);
        MM2=cos(t+k4+2)+cos(t+k4+1+y(dim-1))+cos(t+k4+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        
        Chuk4=mult*multMm1(dim,Chtilde); 
                  
        
        % Boundary Chu(tn+3k10)
        MM1=cos(t+k310+y(1))+cos(t+k310+x(1));
        MM2=cos(t+k310+x(2))+cos(t+k310+y(2))+cos(t+k310);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k310+x(ii));
            MM2=cos(t+k310+x(ii+1))+cos(t+k310+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k310+1+y(1))+cos(t+k310+x(dim));
        MM2=cos(t+k310+1+y(2))+cos(t+k310+1)+cos(t+k310+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;           
        for m=2:dim-1
            MM1=cos(t+k310+y(m));
            MM2=cos(t+k310+y(m+1))+cos(t+k310+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            
            MM1=cos(t+k310+1+y(m));
            MM2=cos(t+k310+1+y(m+1))+cos(t+k310+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end           
        MM1=cos(t+k310+x(1)+1)+cos(t+k310+y(dim));
        MM2=cos(t+k310+x(2)+1)+cos(t+k310+1)+cos(t+k310+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k310+x(ii)+1);
            MM2=cos(t+k310+x(ii+1)+1)+cos(t+k310+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k310+1+y(dim))+cos(t+k310+x(dim)+1);
        MM2=cos(t+k310+2)+cos(t+k310+1+y(dim-1))+cos(t+k310+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        
        Chuk310=mult*multMm1(dim,Chtilde); 
        
        
        % Boundary Chu(tn+3k4)
        MM1=cos(t+k34+y(1))+cos(t+k34+x(1));
        MM2=cos(t+k34+x(2))+cos(t+k34+y(2))+cos(t+k34);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k34+x(ii));
            MM2=cos(t+k34+x(ii+1))+cos(t+k34+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k34+1+y(1))+cos(t+k34+x(dim));
        MM2=cos(t+k34+1+y(2))+cos(t+k34+1)+cos(t+k34+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;           
        for m=2:dim-1
            MM1=cos(t+k34+y(m));
            MM2=cos(t+k34+y(m+1))+cos(t+k34+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            
            MM1=cos(t+k34+1+y(m));
            MM2=cos(t+k34+1+y(m+1))+cos(t+k34+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end           
        MM1=cos(t+k34+x(1)+1)+cos(t+k34+y(dim));
        MM2=cos(t+k34+x(2)+1)+cos(t+k34+1)+cos(t+k34+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k34+x(ii)+1);
            MM2=cos(t+k34+x(ii+1)+1)+cos(t+k34+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k34+1+y(dim))+cos(t+k34+x(dim)+1);
        MM2=cos(t+k34+2)+cos(t+k34+1+y(dim-1))+cos(t+k34+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        
        Chuk34=mult*multMm1(dim,Chtilde);
        
                              
        % Boundary D_h F-ut.
        Dhtilde(1,1)=2*cos(t+y(1))+2*cos(t+x(1));
        for ii=2:dim-1
            Dhtilde(ii,1)=2*cos(t+x(ii));
        end
        Dhtilde(dim,1)=2*cos(t+1+x(1))+2*cos(t+x(dim));        
        for m=2:dim-1
            Dhtilde(1+(m-1)*dim,1)=2*cos(t+y(m));
            Dhtilde(dim+(m-1)*dim,1)=2*cos(t+1+y(m));            
        end        
        Dhtilde(1+(dim-1)*dim,1)=2*cos(t+y(dim))+2*cos(t+x(1)+1);
        for ii=2:dim-1
            Dhtilde(ii+(dim-1)*dim,1)=2*cos(t+x(ii)+1);            
        end
        Dhtilde(dim2,1)=2*cos(t+1+y(dim))+2*cos(t+x(dim)+1); 
        
        Dhfu=multMm1(dim,Dhtilde);                 
       
          
        % Boundary Dh f(t+k/2,u(t+k/2))-ut(t+k2))         
        Dhtilde(1,1)=2*cos(t+k2+y(1))+2*cos(t+k2+x(1));
        for ii=2:dim-1
            Dhtilde(ii,1)=2*cos(t+k2+x(ii));
        end
        Dhtilde(dim,1)=2*cos(t+k2+1+x(1))+2*cos(t+k2+x(dim));        
        for m=2:dim-1
            Dhtilde(1+(m-1)*dim,1)=2*cos(t+k2+y(m));
            Dhtilde(dim+(m-1)*dim,1)=2*cos(t+k2+1+y(m));            
        end        
        Dhtilde(1+(dim-1)*dim,1)=2*cos(t+k2+y(dim))+2*cos(t+k2+x(1)+1);
        for ii=2:dim-1
            Dhtilde(ii+(dim-1)*dim,1)=2*cos(t+k2+x(ii)+1);            
        end
        Dhtilde(dim2,1)=2*cos(t+k2+1+y(dim))+2*cos(t+k2+x(dim)+1);

        Dhfuk2=multMm1(dim,Dhtilde);
             
        % Boundary Dh f(t+k,u(t+k))-ut(t+k)) 
        Dhtilde(1,1)=2*cos(t+k+y(1))+2*cos(t+k+x(1));
        for ii=2:dim-1
            Dhtilde(ii,1)=2*cos(t+k+x(ii));
        end
        Dhtilde(dim,1)=2*cos(t+k+1+x(1))+2*cos(t+k+x(dim));        
        for m=2:dim-1
            Dhtilde(1+(m-1)*dim,1)=2*cos(t+k+y(m));
            Dhtilde(dim+(m-1)*dim,1)=2*cos(t+k+1+y(m));            
        end
        Dhtilde(1+(dim-1)*dim,1)=2*cos(t+k+y(dim))+2*cos(t+k+x(1)+1);
        for ii=2:dim-1
            Dhtilde(ii+(dim-1)*dim,1)=2*cos(t+k+x(ii)+1);            
        end
        Dhtilde(dim2,1)=2*cos(t+k+1+y(dim))+2*cos(t+k+x(dim)+1);

        Dhfuk=multMm1(dim,Dhtilde);
                        
        
        % Boundary Dh f(t+k3,u(t+k3))-ut(t+k3))         
        Dhtilde(1,1)=2*cos(t+k3+y(1))+2*cos(t+k3+x(1));
        for ii=2:dim-1
            Dhtilde(ii,1)=2*cos(t+k3+x(ii));
        end
        Dhtilde(dim,1)=2*cos(t+k3+1+x(1))+2*cos(t+k3+x(dim));        
        for m=2:dim-1
            Dhtilde(1+(m-1)*dim,1)=2*cos(t+k3+y(m));
            Dhtilde(dim+(m-1)*dim,1)=2*cos(t+k3+1+y(m));            
        end        
        Dhtilde(1+(dim-1)*dim,1)=2*cos(t+k3+y(dim))+2*cos(t+k3+x(1)+1);
        for ii=2:dim-1
            Dhtilde(ii+(dim-1)*dim,1)=2*cos(t+k3+x(ii)+1);            
        end
        Dhtilde(dim2,1)=2*cos(t+k3+1+y(dim))+2*cos(t+k3+x(dim)+1);

        Dhfuk3=multMm1(dim,Dhtilde);
        
        % Boundary Dh f(t+k4,u(t+k4))-ut(t+k4)) 
        Dhtilde(1,1)=2*cos(t+k4+y(1))+2*cos(t+k4+x(1));
        for ii=2:dim-1
            Dhtilde(ii,1)=2*cos(t+k4+x(ii));
        end
        Dhtilde(dim,1)=2*cos(t+k4+1+x(1))+2*cos(t+k4+x(dim));        
        for m=2:dim-1
            Dhtilde(1+(m-1)*dim,1)=2*cos(t+k4+y(m));
            Dhtilde(dim+(m-1)*dim,1)=2*cos(t+k4+1+y(m));            
        end
        Dhtilde(1+(dim-1)*dim,1)=2*cos(t+k4+y(dim))+2*cos(t+k4+x(1)+1);
        for ii=2:dim-1
            Dhtilde(ii+(dim-1)*dim,1)=2*cos(t+k4+x(ii)+1);            
        end
        Dhtilde(dim2,1)=2*cos(t+k4+1+y(dim))+2*cos(t+k4+x(dim)+1);

        Dhfuk4=multMm1(dim,Dhtilde);            

                
        % Boundary Dh f(t+k310,u(t+k310))-ut(t+k310))         
        Dhtilde(1,1)=2*cos(t+k310+y(1))+2*cos(t+k310+x(1));
        for ii=2:dim-1
            Dhtilde(ii,1)=2*cos(t+k310+x(ii));
        end
        Dhtilde(dim,1)=2*cos(t+k310+1+x(1))+2*cos(t+k310+x(dim));        
        for m=2:dim-1
            Dhtilde(1+(m-1)*dim,1)=2*cos(t+k310+y(m));
            Dhtilde(dim+(m-1)*dim,1)=2*cos(t+k310+1+y(m));            
        end        
        Dhtilde(1+(dim-1)*dim,1)=2*cos(t+k310+y(dim))+2*cos(t+k310+x(1)+1);
        for ii=2:dim-1
            Dhtilde(ii+(dim-1)*dim,1)=2*cos(t+k310+x(ii)+1);            
        end
        Dhtilde(dim2,1)=2*cos(t+k310+1+y(dim))+2*cos(t+k310+x(dim)+1);

        Dhfuk310=multMm1(dim,Dhtilde);
        
        
        % Boundary Dh f(t+k34,u(t+k34))-ut(t+k34))         
        Dhtilde(1,1)=2*cos(t+k34+y(1))+2*cos(t+k34+x(1));
        for ii=2:dim-1
            Dhtilde(ii,1)=2*cos(t+k34+x(ii));
        end
        Dhtilde(dim,1)=2*cos(t+k34+1+x(1))+2*cos(t+k34+x(dim));        
        for m=2:dim-1
            Dhtilde(1+(m-1)*dim,1)=2*cos(t+k34+y(m));
            Dhtilde(dim+(m-1)*dim,1)=2*cos(t+k34+1+y(m));            
        end        
        Dhtilde(1+(dim-1)*dim,1)=2*cos(t+k34+y(dim))+2*cos(t+k34+x(1)+1);
        for ii=2:dim-1
            Dhtilde(ii+(dim-1)*dim,1)=2*cos(t+k34+x(ii)+1);            
        end
        Dhtilde(dim2,1)=2*cos(t+k34+1+y(dim))+2*cos(t+k34+x(dim)+1);

        Dhfuk34=multMm1(dim,Dhtilde);
        
        
        
        % Evaluations of function h       
        for ii=1:dim
            for jj=1:dim
                pos=ii+(jj-1)*dim;               
                point=t+x(ii)+y(jj);              
                vech(pos,1)=-sin(point)+2*cos(point)-cos(point)^2;
                vechk(pos,1)=-sin(point+k)+2*cos(point+k)-cos(point+k)^2; 
                vechk2(pos,1)=-sin(point+k2)+2*cos(point+k2)-cos(point+k2)^2;
                vechk3(pos,1)=-sin(point+k3)+2*cos(point+k3)-cos(point+k3)^2;  
                vechk4(pos,1)=-sin(point+k4)+2*cos(point+k4)-cos(point+k4)^2; 
                vechk310(pos,1)=-sin(point+k310)+2*cos(point+k310)-cos(point+k310)^2;  
                vechk34(pos,1)=-sin(point+k34)+2*cos(point+k34)-cos(point+k34)^2;                  
            end
        end
                               
        % Stage K2      
        G1=U.^2+vech+Chu+Dhfu;         
        vecb2(:,1)=U;
        vecb2(:,2)=G1;                
        K2=phipmM9points(k2,JJ,vecb2,10^(-13),1,2);          
                      
        % Stages K3 and K4
        D2=K2.^2+vechk2+Chuk2+Dhfuk2-G1;       
        vecb3(:,1)=U;
        vecb3(:,2)=G1;
        vecb3(:,3)=2*D2/k;        
        tie(1)=k3;
        tie(2)=k2;        
        Uaux=phipm_simul_iom9points(tie,JJ,vecb3,10^(-13),1,2);
        K4=Uaux(:,1);
        K3=Uaux(:,2);
                
        
        % Stages K5, K6 and K7
        D3=K3.^2+vechk2+Chuk2+Dhfuk2-G1; 
        D4=K4.^2+vechk3+Chuk3+Dhfuk3-G1;                     
        vecb4(:,1)=U;
        vecb4(:,2)=G1;
        vecb4(:,3)=(-4*D3+9*D4)/k;
        vecb4(:,4)=(24*D3-36*D4)/k^2;          
        tie2(1)=k4;
        tie2(2)=k3;
        tie2(3)=k2;
        Uaux2=phipm_simul_iom9points(tie2,JJ,vecb4,10^(-13),1,2);
        K7=Uaux2(:,1);
        K6=Uaux2(:,2);
        K5=Uaux2(:,3);
                
        
        % Stages K8, K9 and K10         
        D5=K5.^2+vechk2+Chuk2+Dhfuk2-G1;   
        D6=K6.^2+vechk3+Chuk3+Dhfuk3-G1;    
        D7=K7.^2+vechk4+Chuk4+Dhfuk4-G1;          
        vecb5(:,1)=U;
        vecb5(:,2)=G1;
        vecb5(:,3)=(4*D5-27*D6+32*D7)/k;
        vecb5(:,4)=(-56*D5+324*D6-320*D7)/k^2;
        vecb5(:,5)=(288*D5-1296*D6+1152*D7)/k^3;                  
       
        tie2(1)=k310;
        tie2(2)=k34;
        tie2(3)=k;
        Uaux2=phipm_simul_iom9points(tie2,JJ,vecb5,10^(-13),1,2);
        K8=Uaux2(:,1);
        K9=Uaux2(:,2);
        K10=Uaux2(:,3);
              
       
                       
        % Approximation to the exact solution at time t_{n+1}        
        D8=K8.^2+vechk310+Chuk310+Dhfuk310-G1;      
        D9=K9.^2+vechk34+Chuk34+Dhfuk34-G1;
        D10=K10.^2+vechk+Chuk+Dhfuk-G1;        
        vecb5(:,3)=((500/63)*D8-(32/9)*D9+(9/7)*D10)/k;
        vecb5(:,4)=(-(1000/27)*D8+(832/27)*D9-12*D10)/k^2;
        vecb5(:,5)=((4000/63)*D8-(640/9)*D9+(240/7)*D10)/k^3;        
        U=phipmM9points(k,JJ,vecb5,10^(-13),1,2);
        
        % New value of t
        t=t+k;  
        
    end
    
    % CPU time calculus finishes
    telapsed=toc(tstart)
   
    % ERROR1 contains the exact solution at time T 
    for mm=1:dim
        for ii=1:dim
            ERROR1(ii,mm)=cos(t+x(ii)+y(mm))-U(ii+(mm-1)*dim,1);
        end
    end

    % Error in the infinite norm    
    errdef1=norm(max(abs(ERROR1)),inf);
    % The order in the infinite norm is calculated as log2(e0/errdef1), 
    % with e0 the error obtained with k and errdef1 the error obtained with k/2. 
    % When ll=1, as there is not a previous error, it can't  be calculated. 
    % Two consecutive errors are compared
    if ll==1
        errdef1
        e0=errdef1;
    else
        [errdef1 log2(e0/errdef1)]
        e0=errdef1;   
    end
     
    % The new value of n is 2*n, and the new value of k_n=k/2.
    k=k/2;
    n=2*n;      
end
