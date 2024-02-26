% Program that implements the 4th order method for stiff problems called
% expRK4s6 that is constructed in 
% V. T. Luan, Efficient exponential Runge-Kutta methods of high order:
% construction and implementation, BIT Numerical Mathematics,61 (2021)
% 535-560
% 
% avoiding the order reduction with p=3. 
% The problem that it is solved is
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
uxx0=zeros(dim,1);
uxx1=zeros(dim,1);
uy0y=zeros(dim,1);
uy1y=zeros(dim,1);
ux0y=zeros(dim,1);
ux1y=zeros(dim,1);
uyx0=zeros(dim,1);
uyx1=zeros(dim,1);
uxtx0=zeros(dim,1);
uxtx1=zeros(dim,1);
uyt0y=zeros(dim,1);
uyt1y=zeros(dim,1);
uxt0y=zeros(dim,1);
uxt1y=zeros(dim,1);
uytx0=zeros(dim,1);
uytx1=zeros(dim,1);
ERROR1=zeros(dim,dim);
funu=zeros(dim2,1);
vecb5=zeros(dim2,5);
vecb2=zeros(dim2,4);

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

% The program runs for 9 different values of k, from k=1/2, in order to calculate
% the error and the order of the method
for ll=1:9
    
    k=1/n;
    k2=k/2;
    k3=k/3;
    k56=5*k/6;
    
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
                        
        % Boundaries Chut and D_h Au.          
        MM1=-sin(t+y(1))-sin(t+x(1));
        MM2=-sin(t+x(2))-sin(t+y(2))-sin(t); 
        Dhtilde(1,1)=-2*cos(t+y(1))-2*cos(t+x(1));        
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=-sin(t+x(ii));
            MM2=-sin(t+x(ii+1))-sin(t+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
            Dhtilde(ii,1)=-2*cos(t+x(ii));
        end
        MM1=-sin(t+1+y(1))-sin(t+x(dim));
        MM2=-sin(t+1+y(2))-sin(t+1)-sin(t+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;
        Dhtilde(dim,1)=-2*cos(t+1+y(1))-2*cos(t+x(dim));              
        for m=2:dim-1
            MM1=-sin(t+y(m));
            MM2=-sin(t+y(m+1))-sin(t+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(1+(m-1)*dim,1)=-2*cos(t+y(m));            
            MM1=-sin(t+1+y(m));
            MM2=-sin(t+1+y(m+1))-sin(t+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;                        
            Dhtilde(dim+(m-1)*dim,1)=-2*cos(t+1+y(m)); 
        end           
        MM1=-sin(t+x(1)+1)-sin(t+y(dim));
        MM2=-sin(t+x(2)+1)-sin(t+1)-sin(t+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;        
        Dhtilde(1+(dim-1)*dim,1)=-2*cos(t+y(dim))-2*cos(t+x(1)+1); 
        for ii=2:dim-1
            MM1=-sin(t+x(ii)+1);
            MM2=-sin(t+x(ii+1)+1)-sin(t+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(ii+(dim-1)*dim,1)=-2*cos(t+x(ii)+1);  
        end
        MM1=-sin(t+1+y(dim))-sin(t+x(dim)+1);
        MM2=-sin(t+2)-sin(t+1+y(dim-1))-sin(t+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        Dhtilde(dim2,1)=-2*cos(t+1+y(dim))-2*cos(t+x(dim)+1);         
        Chut=mult*multMm1(dim,Chtilde);  
        DhAu=multMm1(dim,Dhtilde);  
        
        % Boundaries DhAut and ChAut
        MM1=2*sin(t+y(1))+2*sin(t+x(1));
        MM2=2*sin(t+x(2))+2*sin(t+y(2))+2*sin(t);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        Dhtilde(1,1)=MM1;
        for ii=2:dim-1
            MM1=2*sin(t+x(ii));
            MM2=2*sin(t+x(ii+1))+2*sin(t+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
            Dhtilde(ii,1)=MM1;
        end
        MM1=2*sin(t+1+y(1))+2*sin(t+x(dim));
        MM2=2*sin(t+1+y(2))+2*sin(t+1)+2*sin(t+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;
        Dhtilde(dim,1)=MM1;           
        for m=2:dim-1
            MM1=2*sin(t+y(m));
            MM2=2*sin(t+y(m+1))+2*sin(t+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(1+(m-1)*dim,1)=MM1;            
            MM1=2*sin(t+1+y(m));
            MM2=2*sin(t+1+y(m+1))+2*sin(t+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(dim+(m-1)*dim,1)=MM1;  
        end           
        MM1=2*sin(t+x(1)+1)+2*sin(t+y(dim));
        MM2=2*sin(t+x(2)+1)+2*sin(t+1)+2*sin(t+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        Dhtilde(1+(dim-1)*dim,1)=MM1;
        for ii=2:dim-1
            MM1=2*sin(t+x(ii)+1);
            MM2=2*sin(t+x(ii+1)+1)+2*sin(t+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(ii+(dim-1)*dim,1)=MM1;
        end
        MM1=2*sin(t+1+y(dim))+2*sin(t+x(dim)+1);
        MM2=2*sin(t+2)+2*sin(t+1+y(dim-1))+2*sin(t+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        Dhtilde(dim2,1)=MM1;        
        ChAut=mult*multMm1(dim,Chtilde);         
        DhAut=multMm1(dim,Dhtilde);  
                       
        % Boundary Chutt
        MM1=-cos(t+y(1))-cos(t+x(1));
        MM2=-cos(t+x(2))-cos(t+y(2))-cos(t);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=-cos(t+x(ii));
            MM2=-cos(t+x(ii+1))-cos(t+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=-cos(t+1+y(1))-cos(t+x(dim));
        MM2=-cos(t+1+y(2))-cos(t+1)-cos(t+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;           
        for m=2:dim-1
            MM1=-cos(t+y(m));
            MM2=-cos(t+y(m+1))-cos(t+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            
            MM1=-cos(t+1+y(m));
            MM2=-cos(t+1+y(m+1))-cos(t+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end           
        MM1=-cos(t+x(1)+1)-cos(t+y(dim));
        MM2=-cos(t+x(2)+1)-cos(t+1)-cos(t+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=-cos(t+x(ii)+1);
            MM2=-cos(t+x(ii+1)+1)-cos(t+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=-cos(t+1+y(dim))-cos(t+x(dim)+1);
        MM2=-cos(t+2)-cos(t+1+y(dim-1))-cos(t+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;        
        Chutt=mult*multMm1(dim,Chtilde); 
                       
        % Boundary Chuttt
        MM1=sin(t+y(1))+sin(t+x(1));
        MM2=sin(t+x(2))+sin(t+y(2))+sin(t);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=sin(t+x(ii));
            MM2=sin(t+x(ii+1))+sin(t+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=sin(t+1+y(1))+sin(t+x(dim));
        MM2=sin(t+1+y(2))+sin(t+1)+sin(t+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;           
        for m=2:dim-1
            MM1=sin(t+y(m));
            MM2=sin(t+y(m+1))+sin(t+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            
            MM1=sin(t+1+y(m));
            MM2=sin(t+1+y(m+1))+sin(t+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end           
        MM1=sin(t+x(1)+1)+sin(t+y(dim));
        MM2=sin(t+x(2)+1)+sin(t+1)+sin(t+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=sin(t+x(ii)+1);
            MM2=sin(t+x(ii+1)+1)+sin(t+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=sin(t+1+y(dim))+sin(t+x(dim)+1);
        MM2=sin(t+2)+sin(t+1+y(dim-1))+sin(t+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;        
        Chuttt=mult*multMm1(dim,Chtilde);         
        
        % Boundary Dh Autt 
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
        DhAutt=multMm1(dim,Dhtilde);
         
        % Evaluations of function h       
        for ii=1:dim
            for jj=1:dim
                pos=ii+(jj-1)*dim;               
                point=t+x(ii)+y(jj);              
                vech(pos,1)=-sin(point)+2*cos(point)-cos(point)^2;
                vechk2(pos,1)=-sin(point+k2)+2*cos(point+k2)-cos(point+k2)^2;
                vechk3(pos,1)=-sin(point+k3)+2*cos(point+k3)-cos(point+k3)^2; 
                vechk56(pos,1)=-sin(point+k56)+2*cos(point+k56)-cos(point+k56)^2;      
            end
        end
                               
        % Stage K2     
        Phfk1=U.^2+vech;
        G1=Phfk1+Chu-DhAu;            
        vecb2(:,1)=U;
        vecb2(:,2)=G1; 
        vecb2(:,3)=Chut-DhAut;
        vecb2(:,4)=ChAut;
        K2=phipmM9points(k2,JJ,vecb2,10^(-13),1,2);      
                      
        % Stages K3 and K4         
        Phfk2=K2.^2+vechk2;         
        vecb2(:,3)=2*(Phfk2-Phfk1)/k+Chut-DhAut; 
        vecb2(:,4)=Chutt;          
        tie(1)=k3;
        tie(2)=k2;        
        Uaux=phipm_simul_iom9points(tie,JJ,vecb2,10^(-13),1,2);
        K4=Uaux(:,1);
        K3=Uaux(:,2);
 
        % Stages K5 and K6
        Phfk3=K3.^2+vechk2;
        Phfk4=K4.^2+vechk3;   
        vecb2(:,3)=(-5*Phfk1-4*Phfk3+9*Phfk4)/k+Chut-DhAut;
        vecb2(:,4)=(12*Phfk1+24*Phfk3-36*Phfk4)/k^2+Chutt;   
        tie(1)=k3;
        tie(2)=k56;        
        Uaux=phipm_simul_iom9points(tie,JJ,vecb2,10^(-13),1,2);
        K6=Uaux(:,1);
        K5=Uaux(:,2);
                     
        
        % Approximation to the exact solution at time t_{n+1} 
        Phfk5=K5.^2+vechk56;
        Phfk6=K6.^2+vechk3;         
        vecb5(:,1)=U;
        vecb5(:,2)=G1;        
        vecb5(:,3)=(-21*Phfk1/5-4*Phfk5/5+5*Phfk6)/k+Chut-DhAut;
        vecb5(:,4)=(36*Phfk1/5+24*Phfk5/5-12*Phfk6)/k^2+Chutt-DhAutt;        
        vecb5(:,5)=Chuttt;             
           
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
