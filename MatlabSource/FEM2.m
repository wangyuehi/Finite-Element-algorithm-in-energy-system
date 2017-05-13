%--------------Written by Behzad Jandaghi, University of Alberta-----------
%-----------------Linear Induction Machine Finite Element Code-------------

clc;
clear all;

%----------------------------initial information---------------------------

NNO=534;    % number of nodes-----Number of nodes changes in each time step due to supplementary nodes but their contribution will be assembled to other nodes so, the matrix size will not be affects by the supplementary nodes
NEL=1008;   % number of elements-----which are fixed in moving structures
ap=0;       % antiperiodicity index for reversing the current direction in each position reset
vel=26.82;  % speed (m/s)
L=0.2667;   % machine depth (m)-----used in force calculation section only

%--------------------------------Mesh--------------------------------------

load mesh; % mesh array that indicates the node numbers of each element-----the moving band related elements will change in each time step
load X;    % array that indicates the X coordinate of the nodes
load Y;    % array that indicates the Y coordinate of the nodes

%--------------------------Time Dependency Modeling------------------------
tic;
dt=1e-3;                       %time step
T=0.001;                           %simulation time
A=zeros(600,round(T/dt+1));    %the solution matrix
A(:,1)=0.001;                  %initial solution
At=zeros(600,1);               % temporary At used to expand the result before saving to A
F_thrust=zeros(round(T/dt+1),1);
F_Normal=zeros(round(T/dt+1),1);

%------------------------meshes speed and conductivity---------------------

%indicates the speed of each element, nonzero for the primary part
v=zeros(NEL,1);
for i=[52:152,220:320,388:488,556:656,724:824,892:992]  %moving part mesh numbers
    v(i)=vel;
end

%indicates the material conductivity of each element, nonzero only for the aluminium sheet 
Sigma=zeros(NEL,1);
for i=[11:26,179:194,347:362,515:530,683:698,851:866]
    Sigma(i)=3.25*10^7;
end

%----------------------------time dependency loop--------------------------

%actually, t+dt:0.001:dt:5, t is the previous time not the current time step
tt=1;                 % a signal for debugging
for t=0:dt:T       %basd on backward euler, we are in t+dt             
        
    
    
    
    %------------------------Moving Band Modeling------------------------------
    
    % changing the coordinates of the moving air gap nodes
    for i=[25:86,111:172,197:258,283:344,369:430,455:516,522:534] % moving nodes of the the moving part
        X(i)=X(i)+vel*dt;
    end
    
    % restarting the position of the rotor when it reaches the end of stator
    if X(25)>0.1524
        ap=ap+1; % the signal for reversing the current direction after each position reset to satisfy the antiperiodicity properties
        for i=[25:86,111:172,197:258,283:344,369:430,455:516,522:534]
            X(i)=X(i)-0.1524;
        end
    end
    
    % finding the appropriate elements in the airgap (excluded the supplementary nodes)
    mbu=[25:32,111:118,197:204,283:290,369:376,455:462,522]; %moving band upper nodes
    mbd=[17:24,103:110,189:196,275:282,361:368,447:454,521]; %moving band downer nodes
    remesh=[153:168,321:336,489:504,657:672,825:840,993:1008]; %moving band meshes
    n_supp=0;
    N_m=49;            % indicates the number of moving band nodes on upper and lower sides
    for i=1:N_m-1      %remeshing the airgap
        if X(25)>=X(mbd(i)) && X(25)<X(mbd(i+1))
            for j=1:N_m-n_supp-1
                mesh(remesh(2*j-1),:)=[mbu(j) mbd(j+n_supp+1) mbd(j+n_supp)];
                mesh(remesh(2*j),:)=[mbu(j) mbu(j+1) mbd(j+n_supp+1)];
            end
            break; % as soon as finding appropriate mesh the for loop should be ended
        else
            n_supp=n_supp+1; % indicates the number of supplementary nodes
        end
    end
    
    % defining the supplementary nodes positions using out nodes due to movement
    for i=N_m-n_supp+1:N_m
        X(i+NNO+n_supp-N_m)=X(mbu(i-1))-0.1524;
        Y(i+NNO+n_supp-N_m)=0.0318;
    end
    
    % defining the new elements related to the supplementary nodes
    if n_supp>0
        
        for i=1:n_supp-1
            mesh(remesh(2*(i+N_m-n_supp-1)-1),:)=[NNO+i mbd(i+1) mbd(i)];
            mesh(remesh(2*(i+N_m-n_supp-1)),:)=[NNO+i NNO+i+1 mbd(i+1)];
        end
        mesh(remesh(2*(N_m-1)-1),:)=[NNO+n_supp mbd(n_supp+1) mbd(n_supp)];
        mesh(remesh(2*(N_m-1)),:)=[NNO+n_supp mbu(1) mbd(n_supp+1)];
    end
    
    %------------------------Specifying the Source-----------------------------
    
    Acu=1.122e-3/2;                      % copper cross section area of a slot divided by 2 for each winding (m2)
    f=92;                                %supply frequency (Hz)
    
    Amp=sqrt(2)*3*1541/(4*Acu);          %supply current density amplitude (A)  J=NI/(c*A) where c is the number of parallel circuits
    
    Ja=Amp*sin(2*pi*f*(t+dt));           %Current density of phase a
    Jb=Amp*sin(2*pi*f*(t+dt)-4*pi/3);
    Jc=Amp*sin(2*pi*f*(t+dt)-2*pi/3);
    
    J=zeros(NEL,1);                      %Current density, only available in primary mover windings
    
    % for considering the reverse current when the mover reaches the end and reset to the begin
    
    for i=[114:136,137:152,305:320]      % element numbers where Ia exists
        J(i)=(-1)^ap*Ja;                 % first slot up down+second slot up
    end
    
    for i=954:976
        J(i)=-1*(-1)^ap*Ja;
    end
    
    for i=[618:640,786:808,809:824,977:992]
        J(i)=(-1)^ap*Jb;
    end
    
    for i=[282:304,450:472,473:488,641:656]
        J(i)=-1*(-1)^ap*Jc;
    end
    
    %--------------------------------------------------------------------------
    %N_supp is now known, so the dimension of the matrices are known
    
    An(:,1)=ones(NNO+n_supp,1)*0.001;       % the first iteration, initial solution of the corresponding time step
    
    for I=1:20                %magnetic material nonlinearity loop, I is the iteration number
        
        %------------------------Specifying the Materials--------------------------
        
        %indicates the material permieability of each element, 1 for air, copper and the aluminium sheet regions
        
        nu=ones(NEL,1)/(4*pi*1e-7);
        dnoodB2=zeros(NEL,1);
        
        for i=[1:10,169:178,337:346,505:514,673:682,841:850]                   % the secondary core element numbers 1010 steel
            N1=mesh(i,1);                                                      % Global node numbers of that element
            N2=mesh(i,2);
            N3=mesh(i,3);
            DET=X(N2)*Y(N3)+X(N1)*Y(N2)+X(N3)*Y(N1)-X(N1)*Y(N3)-X(N3)*Y(N2)-X(N2)*Y(N1);
            %B=curl(A)=dA/dy i - dA/dx j , where A=a1+a2x+a3y, so dA/dy= a3 and dA/dx=a2
            a2=(An(N2,I)*Y(N3)-An(N3,I)*Y(N2)-An(N1,I)*(Y(N3)-Y(N2))+Y(N1)*(An(N3,I)-An(N2,I)))/DET;
            a3=(An(N3,I)*X(N2)-An(N2,I)*X(N3)-X(N1)*(An(N3,I)-An(N2,I))+An(N1,I)*(X(N3)-X(N2)))/DET;
            B(i)=sqrt(a2^2+a3^2);
            %preventing B from having zero value
            if B(i)==0
                B(i)=1e-12;
            end
            %magnetic saturation curve "1010 steel"
            if B(i)<0.6
                nu(i)=135.28;
                dnoodB2(i)=0;
            else
                nu(i)=79.57*((B(i)-0.6)^2+1.7);
                dnoodB2(i)=79.57;
            end
        end
        
        for i=[52:99,220:267,388:435,556:603,724:771,892:939]                     % the primary core element numbers
            N1=mesh(i,1);                                                       % Global node number of that element
            N2=mesh(i,2);
            N3=mesh(i,3);
            DET=X(N2)*Y(N3)+X(N1)*Y(N2)+X(N3)*Y(N1)-X(N1)*Y(N3)-X(N3)*Y(N2)-X(N2)*Y(N1);
            %B=curl(A)=dA/dy i - dA/dx j , where A=a1+a2x+a3y, so dA/dy= a3 and dA/dx=a2
            a2=(An(N2,I)*Y(N3)-An(N3,I)*Y(N2)-An(N1,I)*(Y(N3)-Y(N2))+Y(N1)*(An(N3,I)-An(N2,I)))/DET;
            a3=(An(N3,I)*X(N2)-An(N2,I)*X(N3)-X(N1)*(An(N3,I)-An(N2,I))+An(N1,I)*(X(N3)-X(N2)))/DET;
            B(i)=sqrt(a2^2+a3^2);
            %preventing B from having zero value
            if B(i)==0
                B(i)=1e-12;
            end
            %magnetic saturation curve "M19"
            if B(i)<0.6
                nu(i)=135.28;
                dnoodB2(i)=0;
            else
                nu(i)=79.57*((B(i)-0.6)^2+1.7);
                dnoodB2(i)=79.57;
            end
        end
        
        %-----------------------------Matrix equation------------------------------
        %matrices should be reset to zero for each iteration
        SS=zeros(NNO+n_supp,NNO+n_supp); %stifness matrix
        SJ=zeros(NNO+n_supp,NNO+n_supp); % Jakoobian Matrix
        b=zeros(NNO+n_supp,1);    %source term
        
        %-----------------------stiffness matrix formation-------------------------
        
        for i=1:NEL
            
            N1=mesh(i,1); % Global node number of that element
            N2=mesh(i,2);
            N3=mesh(i,3);
            
            Q1=Y(N2)-Y(N3);
            Q2=Y(N3)-Y(N1);
            Q3=Y(N1)-Y(N2);
            R1=X(N3)-X(N2);
            R2=X(N1)-X(N3);
            R3=X(N2)-X(N1);
            
            DET=X(N2)*Y(N3)+X(N1)*Y(N2)+X(N3)*Y(N1)-X(N1)*Y(N3)-X(N3)*Y(N2)-X(N2)*Y(N1);
            
            COEFF1=nu(i)/(DET*2);
            COEFF2=Sigma(i)*DET/(12*dt);
            COEFF3=Sigma(i)*v(i)/6;
            COEFF4=J(i)*DET/6;
            
            %----stiffness matrix for each element corresponding to the main term------
            
            S1(1,1)=COEFF1*(Q1*Q1+R1*R1);
            S1(1,2)=COEFF1*(Q1*Q2+R1*R2);
            S1(1,3)=COEFF1*(Q1*Q3+R1*R3);
            S1(2,1)=S1(1,2);
            S1(2,2)=COEFF1*(Q2*Q2+R2*R2);
            S1(2,3)=COEFF1*(Q2*Q3+R2*R3);
            S1(3,1)=S1(1,3);
            S1(3,2)=S1(2,3);
            S1(3,3)=COEFF1*(Q3*Q3+R3*R3);
            
            for k=1:3
                for n=1:3
                    J1(k,n)=S1(k,n)+4/(nu(i)^2*DET)*(dnoodB2(i))*sum(S1(k,1)*An(N1,I)+S1(k,2)*An(N2,I)+S1(k,3)*An(N3,I))*sum(S1(n,1)*An(N1,I)+S1(n,2)*An(N2,I)+S1(n,3)*An(N3,I));            %where is A1,A2,A3????????????
                end
            end
            
            %--stiffness matrix for each element corresponding to the conducting part--
            
            S2(1,1)=COEFF2;
            S2(1,2)=COEFF2*0.5;
            S2(1,3)=COEFF2*0.5;
            S2(2,1)=S2(1,2);
            S2(2,2)=COEFF2;
            S2(2,3)=COEFF2*0.5;
            S2(3,1)=S2(1,3);
            S2(3,2)=S2(2,3);
            S2(3,3)=COEFF2;
            
            %----stiffness matrix for each element corresponding to the moving part----
            
            S3(1,1)=COEFF3*Q1;
            S3(1,2)=COEFF3*Q2;
            S3(1,3)=COEFF3*Q3;
            S3(2,1)=S3(1,1);
            S3(2,2)=S3(1,2);
            S3(2,3)=S3(1,3);
            S3(3,1)=S3(1,1);
            S3(3,2)=S3(1,2);
            S3(3,3)=S3(1,3);
            
            %-------source term corresponding to the eddy current current source-------
            
            b1(1)=COEFF2*(A(N1,round(t/dt+1))+0.5*A(N2,round(t/dt+1))+0.5*A(N3,round(t/dt+1)));
            b1(2)=COEFF2*(0.5*A(N1,round(t/dt+1))+A(N2,round(t/dt+1))+0.5*A(N3,round(t/dt+1)));
            b1(3)=COEFF2*(0.5*A(N1,round(t/dt+1))+0.5*A(N2,round(t/dt+1))+A(N3,round(t/dt+1)));
            
            %---------source term corresponding to the external current source---------
            
            b2(1)=COEFF4;
            b2(2)=COEFF4;
            b2(3)=COEFF4;
            
            %     %#################################
            %     if tt==6
            %         stop
            %     end
            %     %#################################
            
            %----------------construction of the global stiffness matrix---------------
            
            for i1=1:3
                b(mesh(i,i1))=b(mesh(i,i1))+b1(i1)+b2(i1); % source term
                for i2=1:3
                    SS(mesh(i,i1),mesh(i,i2))=SS(mesh(i,i1),mesh(i,i2))+S1(i1,i2)+S2(i1,i2)+S3(i1,i2); %inserting local stiffness matrix to global stiffness matrix
                end
            end
            
            %----------------construction of the global Jakoonian matrix---------------
            %the local jakoobian matrices J2 and J3 are equal to S2 and S3
            %SJ(mesh(i,i1),mesh(i,i2))=0;
            for i1=1:3
                for i2=1:3
                    SJ(mesh(i,i1),mesh(i,i2))=SJ(mesh(i,i1),mesh(i,i2))+J1(i1,i2)+S2(i1,i2)+S3(i1,i2); %inserting local stiffness matrix to global Jakoobian matrix
                end
            end
            
        end
        
        %--------------------------imposing boundary conditions--------------------
        
        % anti periodicity boundary condition on the right hand side of the machine
        
        L_P=[1,3,5,13,17,25,33,43,48,55,59,63,67,71,75,79,83,85,NNO+1:NNO+n_supp]; % node number on the anitiperiodicity boundary left
        R_P=[517:534,mbu(N_m-n_supp:N_m-1)]; % node number on the anitiperiodicity boundary right
        
        % for stifness matrix
        
        for i=1:18+n_supp
            for j=1:18+n_supp
                if i==j
                    SS(R_P(i),R_P(j))=SS(R_P(i),R_P(j))+SS(L_P(i),L_P(j));
                else
                    SS(R_P(i),R_P(j))=SS(R_P(i),R_P(j))+SS(L_P(i),L_P(j));
                end
            end
        end
        
        for i=1:18+n_supp
            b(R_P(i))=b(R_P(i))-b(L_P(i));
            for j=1:NNO+n_supp
                if ismember(j,L_P)==0
                    SS(R_P(i),j)=SS(R_P(i),j)-SS(L_P(i),j);
                    SS(j,R_P(i))=SS(j,R_P(i))-SS(j,L_P(i));
                end
            end
        end
        
        % for Jakoobian matrix
        
        for i=1:18+n_supp
            for j=1:18+n_supp
                if i==j
                    SJ(R_P(i),R_P(j))=SJ(R_P(i),R_P(j))+SJ(L_P(i),L_P(j));
                else
                    SJ(R_P(i),R_P(j))=SJ(R_P(i),R_P(j))+SJ(L_P(i),L_P(j));
                end
            end
        end
        
        for i=1:18+n_supp
            for j=1:NNO+n_supp
                if ismember(j,L_P)==0
                    SJ(R_P(i),j)=SJ(R_P(i),j)-SJ(L_P(i),j);
                    SJ(j,R_P(i))=SJ(j,R_P(i))-SJ(j,L_P(i));
                end
            end
        end
        
        %BEL=[1,2,6,5,12,11,20,28,27,remesh(96-2*n_supp+1):remesh(96),remesh(1),remesh(2),40,39,52,53,54,60,61,62,63,68,69,70,74,75,76,81,82,83,84,89,90,91,92,97,96] % element numbers on the boundary that have at least one node on the boundary
        
        % Drichlet boundary condition (A=0) on top and bottom sides of the geometry
        % in order to make BC nodes to have zero potential we have to make zero the corresponding nodes potential of 1- An 2- SJ\(SS*An(:,I)-b) which results in:
        for i=[1,2,87,88,173,174,259,260,345,346,431,432,517,85,86,171,172,257,258,343,344,429,430,515,516,534] %node numbers at which A=0
            SS(i,:)=0;
            SS(i,i)=1;
            SJ(i,:)=0;
            SJ(i,i)=1;
            b(i)=0;
        end
        
        % up to this point the size of matrices are NNO+n_supp
        % removing the extra rows and columns of the stifness and jakoobian matrices related to the anitiperiodicity boundary nodes of the left side,supplementary nodes
        
        k1=0;
        for i=[1,3,5,13,17,25,33,43,48,55,59,63,67,71,75,79,83,85,NNO+1:NNO+n_supp]
            SS(i-k1,:)=[];
            SS(:,i-k1)=[];
            b(i-k1)=[];
            SJ(i-k1,:)=[];
            SJ(:,i-k1)=[];
            An(i-k1,:)=[];
            k1=k1+1;
        end
        
        %------------------------------NR solution---------------------------------
        % SS,SJ and b are now 516*516
        
        % Naive solver code
        num=SJ;
        den=(SS*An(:,I)-b);
        frac=num\den;
        ns=516;
        
        num(:,ns+1)=den;
        for is=1:ns
            %Step 2: make diagonal elements into 1.0
            num(is,is+1:ns+1) = num(is,is+1:ns+1)/num(is,is);
            %Step 3: make all elements below diagonal into 0
            for js=is+1:ns
                num(js,is+1:ns+1) = num(js,is+1:ns+1) - num(js,is)*num(is,is+1:ns+1);
            end
        end
        %Step 4: begin back substitution
        for js=ns-1:-1:1
            num(js,ns+1) = num(js,ns+1) - num(js,js+1:ns)*num(js+1:ns,ns+1);
        end
        %return solution
        frac=num(:,ns+1);
        
        An(:,I+1)=An(:,I)-frac;
        
        % Convergency criteria
        
        if norm(An(:,I+1)-An(:,I))/norm(An(:,I))<5e-3    %convergency condition
             display(I);
            break % when the solution is found, the for loop will break with a 516*1 An
        end
        
        %###################################### considered the last solution as the correct solution without satisfying the convergency criteria
        if I==20 %limit of NR iterations
            %Disp('Not converged')
            break
        end
        %######################################
        
        % expanding the An to the NNO+N_supp matrix required for the next iteration
        k2=1;
        for i=[1,3,5,13,17,25,33,43,48,55,59,63,67,71,75,79,83,85,NNO+1:NNO+n_supp]
            An(i+1:516+k2,:)=An(i:515+k2,:);
            An(i,:)=0;
            k2=k2+1;
        end
        % these two loops should be seperate as it is, should not be combined
        k3=1;
        for i=[1,3,5,13,17,25,33,43,48,55,59,63,67,71,75,79,83,85,NNO+1:NNO+n_supp]
            An(i,:)=-1*An(R_P(k3),:);
            k3=k3+1;
        end
        
%         display (I)
        
    end   %end of NR loop
    
    At(1:516,1)=An(:,I+1);
    
    % expanding the At to the NNO+N_supp matrix required for the next time step
    
    k2=1;
    for i=[1,3,5,13,17,25,33,43,48,55,59,63,67,71,75,79,83,85,NNO+1:NNO+n_supp]
        At(i+1:516+k2,:)=At(i:515+k2,:);
        At(i,:)=0;
        k2=k2+1;
    end
    
    k3=1;
    for i=[1,3,5,13,17,25,33,43,48,55,59,63,67,71,75,79,83,85,NNO+1:NNO+n_supp]
        At(i,:)=-1*At(R_P(k3),:);
        k3=k3+1;
    end
    
    A(1:NNO+n_supp,round(t/dt+2))=At(1:NNO+n_supp,1);% saving expanded At to A as the final solution of each time step
    
    %----------------------Calculation of Forces-------------------------------
    ii=[27,28,28,29,30,31,31,32,33,34,34,35,36,37,37,38];
    
    for i=[ii,ii+168,ii+336,ii+504,ii+672,ii+840]                   % the element numbers corresponding to the force calculation band
        N1=mesh(i,1);                                                      % Global node number of that element
        N2=mesh(i,2);
        N3=mesh(i,3);
        DET=X(N2)*Y(N3)+X(N1)*Y(N2)+X(N3)*Y(N1)-X(N1)*Y(N3)-X(N3)*Y(N2)-X(N2)*Y(N1);
        %B=curl(A)=dA/dy i - dA/dx j , where A=a1+a2x+a3y, so dA/dy= a3 and dA/dx=a2
        a2=(A(N2,round(t/dt+2))*Y(N3)-A(N3,round(t/dt+2))*Y(N2)-A(N1,round(t/dt+2))*(Y(N3)-Y(N2))+Y(N1)*(A(N3,round(t/dt+2))-A(N2,round(t/dt+2))))/DET;
        a3=(A(N3,round(t/dt+2))*X(N2)-A(N2,round(t/dt+2))*X(N3)-X(N1)*(A(N3,round(t/dt+2))-A(N2,round(t/dt+2)))+A(N1,round(t/dt+2))*(X(N3)-X(N2)))/DET;
        F_thrust(round(t/dt+2),1)=F_thrust(round(t/dt+2))+21*L*a2*a3/(4*pi*1e-7)*0.0254/16;
        F_Normal(round(t/dt+2),1)=F_Normal(round(t/dt+2))+21*L*(a2^2-a3^2)/(2*4*pi*1e-7)*0.0254/16;
    end
    
     display(t);
    
    %#################################
    %     tt=tt+1
    %     if tt==70
    %         break
    %     end
    %#################################
    
    clear An
end

toc





%{ 

     this is the end of time dependency loop  

%}

% % % % %-------------testing the magnetic field on the secondary core------------- 
% % % % 
% % % % % magnetic field density in the primary mover yoke in respect to position
% % % % k5=1;
% % % % for i=[5:10,173:178,341:346,509:514,677:682,845:850]                     % the secondary core element numbers
% % % %     N1=mesh(i,1);                                                       % Global node number of that element
% % % %     N2=mesh(i,2);
% % % %     N3=mesh(i,3);
% % % %     DET=X(N2)*Y(N3)+X(N1)*Y(N2)+X(N3)*Y(N1)-X(N1)*Y(N3)-X(N3)*Y(N2)-X(N2)*Y(N1);
% % % %     %B=curl(A)=dA/dy i - dA/dx j , where A=a1+a2x+a3y, so dA/dy= a3 and dA/dx=a2
% % % %     a2=(A(N2,2)*Y(N3)-A(N3,2)*Y(N2)-A(N1,2)*(Y(N3)-Y(N2))+Y(N1)*(A(N3,2)-A(N2,2)))/DET;
% % % %     a3=(A(N3,2)*X(N2)-A(N2,2)*X(N3)-X(N1)*(A(N3,2)-A(N2,2))+A(N1,2)*(X(N3)-X(N2)))/DET;
% % % %     XX(k5)=sqrt(a2^2+a3^2);
% % % %     k5=k5+1;
% % % % end
% % % % plot(1:36,XX(1,:))
% % % % 
% % % % % magnetic field density in the secondary yoke in respect to position
% % % % for beh=1:100
% % % % k5=1;
% % % % for i=[5:10,173:178,341:346,509:514,677:682,845:850]                     % the secondary core element numbers
% % % %     N1=mesh(i,1);                                                       % Global node number of that element
% % % %     N2=mesh(i,2);
% % % %     N3=mesh(i,3);
% % % %     DET=X(N2)*Y(N3)+X(N1)*Y(N2)+X(N3)*Y(N1)-X(N1)*Y(N3)-X(N3)*Y(N2)-X(N2)*Y(N1);
% % % %     %B=curl(A)=dA/dy i - dA/dx j , where A=a1+a2x+a3y, so dA/dy= a3 and dA/dx=a2
% % % %     a2=(A(N2,beh)*Y(N3)-A(N3,beh)*Y(N2)-A(N1,beh)*(Y(N3)-Y(N2))+Y(N1)*(A(N3,beh)-A(N2,beh)))/DET;
% % % %     a3=(A(N3,beh)*X(N2)-A(N2,beh)*X(N3)-X(N1)*(A(N3,beh)-A(N2,beh))+A(N1,beh)*(X(N3)-X(N2)))/DET;
% % % %     XX(k5)=sqrt(a2^2+a3^2);
% % % %     k5=k5+1;
% % % % end
% % % % plot(1:36,XX(1,:))
% % % % hold on
% % % % end
% % % % 
% % % % % magnetic flux density of a mesh in respect to time
% % % % k5=1;
% % % % for beh=1:100
% % % %     i=5;                     % the secondary core element numbers
% % % %     N1=mesh(i,1);                                                       % Global node number of that element
% % % %     N2=mesh(i,2);
% % % %     N3=mesh(i,3);
% % % %     DET=X(N2)*Y(N3)+X(N1)*Y(N2)+X(N3)*Y(N1)-X(N1)*Y(N3)-X(N3)*Y(N2)-X(N2)*Y(N1);
% % % %     %B=curl(A)=dA/dy i - dA/dx j , where A=a1+a2x+a3y, so dA/dy= a3 and dA/dx=a2
% % % %     a2=(A(N2,beh)*Y(N3)-A(N3,beh)*Y(N2)-A(N1,beh)*(Y(N3)-Y(N2))+Y(N1)*(A(N3,beh)-A(N2,beh)))/DET;
% % % %     a3=(A(N3,beh)*X(N2)-A(N2,beh)*X(N3)-X(N1)*(A(N3,beh)-A(N2,beh))+A(N1,beh)*(X(N3)-X(N2)))/DET;
% % % %     XX(k5)=sqrt(a2^2+a3^2);
% % % %     k5=k5+1;
% % % % end
% % % % plot(1:100,XX(1,:))
% % % % 
% % % % % magnetic field density in the primary mover yoke in respect to position
% % % % for beh=1:100
% % % % k5=1;
% % % % for i=[96:99,264:267,432:435,600:603,768:771,936:939]                     % the secondary core element numbers
% % % %     N1=mesh(i,1);                                                       % Global node number of that element
% % % %     N2=mesh(i,2);
% % % %     N3=mesh(i,3);
% % % %     DET=X(N2)*Y(N3)+X(N1)*Y(N2)+X(N3)*Y(N1)-X(N1)*Y(N3)-X(N3)*Y(N2)-X(N2)*Y(N1);
% % % %     %B=curl(A)=dA/dy i - dA/dx j , where A=a1+a2x+a3y, so dA/dy= a3 and dA/dx=a2
% % % %     a2=(A(N2,beh)*Y(N3)-A(N3,beh)*Y(N2)-A(N1,beh)*(Y(N3)-Y(N2))+Y(N1)*(A(N3,beh)-A(N2,beh)))/DET;
% % % %     a3=(A(N3,beh)*X(N2)-A(N2,beh)*X(N3)-X(N1)*(A(N3,beh)-A(N2,beh))+A(N1,beh)*(X(N3)-X(N2)))/DET;
% % % %     XX(k5)=sqrt(a2^2+a3^2);
% % % %     k5=k5+1;
% % % % end
% % % % plot(1:24,XX(1,:))
% % % % hold on
% % % % end

