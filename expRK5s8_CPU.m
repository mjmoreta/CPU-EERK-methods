% Program that implements the 5th order method for siff problems and
% vanishing boundary conditions that is constructed in 
% V. T. Luan and A. Ostermann, Explicit exponential Runge-Kutta methods of 
% high order for parabolic problems, J. Comput. Appl. Math. 262 (2014)
% 361-372
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
vechk23=zeros(dim2,1); 
ERROR1=zeros(dim,dim);
funu=zeros(dim2,1);

x=zeros(dim,1);
y=zeros(dim,1);
x=[h:h:1]';
y=x;

% Matrix A_{h,0} is given by M^{-1}A. M^{-1}, is not computed, when
% necessary, a system is solved or phipmM9points is used. Matrices A and M 
% are not computed because they are very large. We use their values different 
% from 0 when necessary

% n is such that the time step size is k=1/n.
n=2;
k=1/n;
Tf=n*k;
h2=h^2;
mult=12/h2;
mult2=h2/12;

% The program runs for 9 different values of k, from k=1/2, in order to calculate
% the error and the order of the method
for ll=1:9
    
    k=1/n;
    k2=k/2;
    k4=k/4;
    k5=k/5;
    k23=2*k/3;
    
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
    for kk=1:1
        
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
                
        % Boundary Chu(tn+k5)
        MM1=cos(t+k5+y(1))+cos(t+k5+x(1));
        MM2=cos(t+k5+x(2))+cos(t+k5+y(2))+cos(t+k5);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k5+x(ii));
            MM2=cos(t+k5+x(ii+1))+cos(t+k5+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k5+1+y(1))+cos(t+k5+x(dim));
        MM2=cos(t+k5+1+y(2))+cos(t+k5+1)+cos(t+k5+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;           
        for m=2:dim-1
            MM1=cos(t+k5+y(m));
            MM2=cos(t+k5+y(m+1))+cos(t+k5+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            
            MM1=cos(t+k5+1+y(m));
            MM2=cos(t+k5+1+y(m+1))+cos(t+k5+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end           
        MM1=cos(t+k5+x(1)+1)+cos(t+k5+y(dim));
        MM2=cos(t+k5+x(2)+1)+cos(t+k5+1)+cos(t+k5+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k5+x(ii)+1);
            MM2=cos(t+k5+x(ii+1)+1)+cos(t+k5+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k5+1+y(dim))+cos(t+k5+x(dim)+1);
        MM2=cos(t+k5+2)+cos(t+k5+1+y(dim-1))+cos(t+k5+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        
        Chuk5=mult*multMm1(dim,Chtilde); 
        
        
        % Boundary Chu(tn+2k3)
        MM1=cos(t+k23+y(1))+cos(t+k23+x(1));
        MM2=cos(t+k23+x(2))+cos(t+k23+y(2))+cos(t+k23);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k23+x(ii));
            MM2=cos(t+k23+x(ii+1))+cos(t+k23+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k23+1+y(1))+cos(t+k23+x(dim));
        MM2=cos(t+k23+1+y(2))+cos(t+k23+1)+cos(t+k23+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;           
        for m=2:dim-1
            MM1=cos(t+k23+y(m));
            MM2=cos(t+k23+y(m+1))+cos(t+k23+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            
            MM1=cos(t+k23+1+y(m));
            MM2=cos(t+k23+1+y(m+1))+cos(t+k23+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end           
        MM1=cos(t+k23+x(1)+1)+cos(t+k23+y(dim));
        MM2=cos(t+k23+x(2)+1)+cos(t+k23+1)+cos(t+k23+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k23+x(ii)+1);
            MM2=cos(t+k23+x(ii+1)+1)+cos(t+k23+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k23+1+y(dim))+cos(t+k23+x(dim)+1);
        MM2=cos(t+k23+2)+cos(t+k23+1+y(dim-1))+cos(t+k23+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        
        Chuk23=mult*multMm1(dim,Chtilde);                       
                      
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
                
        % Boundary Dh f(t+k5,u(t+k5))-ut(t+k5))         
        Dhtilde(1,1)=2*cos(t+k5+y(1))+2*cos(t+k5+x(1));
        for ii=2:dim-1
            Dhtilde(ii,1)=2*cos(t+k5+x(ii));
        end
        Dhtilde(dim,1)=2*cos(t+k5+1+x(1))+2*cos(t+k5+x(dim));        
        for m=2:dim-1
            Dhtilde(1+(m-1)*dim,1)=2*cos(t+k5+y(m));
            Dhtilde(dim+(m-1)*dim,1)=2*cos(t+k5+1+y(m));            
        end        
        Dhtilde(1+(dim-1)*dim,1)=2*cos(t+k5+y(dim))+2*cos(t+k5+x(1)+1);
        for ii=2:dim-1
            Dhtilde(ii+(dim-1)*dim,1)=2*cos(t+k5+x(ii)+1);            
        end
        Dhtilde(dim2,1)=2*cos(t+k5+1+y(dim))+2*cos(t+k5+x(dim)+1);

        Dhfuk5=multMm1(dim,Dhtilde);
                
        % Boundary Dh f(t+k23,u(t+k23))-ut(t+k23))         
        Dhtilde(1,1)=2*cos(t+k23+y(1))+2*cos(t+k23+x(1));
        for ii=2:dim-1
            Dhtilde(ii,1)=2*cos(t+k23+x(ii));
        end
        Dhtilde(dim,1)=2*cos(t+k23+1+x(1))+2*cos(t+k23+x(dim));        
        for m=2:dim-1
            Dhtilde(1+(m-1)*dim,1)=2*cos(t+k23+y(m));
            Dhtilde(dim+(m-1)*dim,1)=2*cos(t+k23+1+y(m));            
        end        
        Dhtilde(1+(dim-1)*dim,1)=2*cos(t+k23+y(dim))+2*cos(t+k23+x(1)+1);
        for ii=2:dim-1
            Dhtilde(ii+(dim-1)*dim,1)=2*cos(t+k23+x(ii)+1);            
        end
        Dhtilde(dim2,1)=2*cos(t+k23+1+y(dim))+2*cos(t+k23+x(dim)+1);

        Dhfuk23=multMm1(dim,Dhtilde);
        
        % Evaluations of function h       
        for ii=1:dim
            for jj=1:dim
                pos=ii+(jj-1)*dim;               
                point=t+x(ii)+y(jj);              
                vech(pos,1)=-sin(point)+2*cos(point)-cos(point)^2;
                vechk2(pos,1)=-sin(point+k2)+2*cos(point+k2)-cos(point+k2)^2;
                vechk(pos,1)=-sin(point+k)+2*cos(point+k)-cos(point+k)^2; 
                vechk4(pos,1)=-sin(point+k4)+2*cos(point+k4)-cos(point+k4)^2; 
                vechk5(pos,1)=-sin(point+k5)+2*cos(point+k5)-cos(point+k5)^2;  
                vechk23(pos,1)=-sin(point+k23)+2*cos(point+k23)-cos(point+k23)^2;  
                
                
            end
        end
                               
        % Stage K2      
        G1=U.^2+vech+Chu+Dhfu;               
        vecb2=zeros(dim2,2);
        vecb2(:,1)=U;
        vecb2(:,2)=G1;                
        K2=phipmM9points(k2,JJ,vecb2,10^(-13),1,1);      
                      
        % Stage K3 
        G2=K2.^2+vechk2+Chuk2+Dhfuk2-G1;        
        vecb3=zeros(dim2,3);
        vecb3(:,1)=U;
        vecb3(:,2)=G1;
        vecb3(:,3)=2*G2/k;        
        K3=phipmM9points(k2,JJ,vecb3,10^(-13),1,1);          
           
        
        % Stage K44 
        G3=K3.^2+vechk2+Chuk2+Dhfuk2-G1;         
        vecb3(:,1)=U;
        vecb3(:,2)=G1;
        vecb3(:,3)=2*G3/k;                     
        K4=phipmM9points(k4,JJ,vecb3,10^(-13),1,1); 
        
        % Stage K5
        G4=K4.^2+vechk4+Chuk4+Dhfuk4-G1;         
        vecb4=zeros(dim2,4);       
        vecb4(:,1)=U;
        vecb4(:,2)=G1;
        vecb4(:,3)=(-2*G3+8*G4)/k;
        vecb4(:,4)=(16*G3-32*G4)/k^2;          
        K5=phipmM9points(k2,JJ,vecb4,10^(-13),1,1);  
  
        % Stage K6                      
        G5=K5.^2+vechk2+Chuk2+Dhfuk2-G1;        
        vecb4(:,1)=U;
        vecb4(:,2)=G1;
        vecb4(:,3)=(8*G4-2*G5)/k;
        vecb4(:,4)=(-32*G4+16*G5)/k^2;           
        K6=phipmM9points(k5,JJ,vecb4,10^(-13),1,1);
               
        % Stage K7                      
        G6=K6.^2+vechk5+Chuk5+Dhfuk5-G1;        
        vecb4(:,1)=U;
        vecb4(:,2)=G1;
        vecb4(:,3)=((-4/3)*G5+(25/3)*G6)/k;
        vecb4(:,4)=((40/3)*G5-(100/3)*G6)/k^2;        
        vecb42=zeros(dim2,4);
        vecb42(:,3)=25*((-20/81)*G4+(5/243)*G5+(125/486)*G6)/k;
        vecb42(:,4)=125*((16/81)*G4-(4/243)*G5-(50/243)*G6)/k^2;       
        K7=phipmM9points(k23,JJ,vecb4,10^(-13),1,1)+phipmM9points(k5,JJ,vecb42,10^(-13),1,1);  
        
        % Stage K8                      
        G7=K7.^2+vechk23+Chuk23+Dhfuk23-G1;        
        vecb5=zeros(dim2,5);
        vecb5(:,1)=U;
        vecb5(:,2)=G1;
        vecb5(:,3)=((-16/3)*G5+(250/21)*G6+(27/14)*G7)/k;
        vecb5(:,4)=((208/3)*G5-(250/3)*G6-27*G7)/k^2;
        vecb5(:,5)=(-240*G5+(1500/7)*G6+(810/7)*G7)/k^3;           
        vecb52=zeros(dim2,5);
        vecb52(:,3)=25*((-4/7)*G5+(25/49)*G6+(27/98)*G7)/k;
        vecb52(:,4)=125*((8/5)*G5-(10/7)*G6-(27/35)*G7)/k^2;
        vecb52(:,5)=625*(-(48/35)*G5+(60/49)*G6+(162/245)*G7)/k^3;        
        vecb53=zeros(dim2,5);
        vecb53(:,3)=(9/4)*((-288/35)*G5+(360/49)*G6+(972/245)*G7)/k;
        vecb53(:,4)=(27/8)*((384/5)*G5-(480/7)*G6-(1296/35)*G7)/k^2;
        vecb53(:,5)=(81/16)*((-1536/7)*G5+(9600/49)*G6+(5184/49)*G7)/k^3;        
        sum1=phipmM9points(k,JJ,vecb5,10^(-13),1,1);        
        sum2=phipmM9points(k5,JJ,vecb52,10^(-13),1,1);
        sum3=phipmM9points(k23,JJ,vecb53,10^(-13),1,1);        
        K8=sum1+sum2+sum3; 
        
        % % Approximation to the exact solution at time t_{n+1}        
        G8=K8.^2+vechk+Chuk+Dhfuk-G1;        
        vecb5(:,3)=((125/14)*G6-(27/14)*G7+G8/2)/k;
        vecb5(:,4)=((-625/14)*G6+(162/7)*G7-(13/2)*G8)/k^2;
        vecb5(:,5)=((1125/14)*G6-(405/7)*G7+(45/2)*G8)/k^3;
        
        U=phipmM9points(k,JJ,vecb5,10^(-13),1,1);
        
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
