function out  = GWO_1D_MPPT_mf(ACT,P)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% variables initialization %%%%%%%%%%%%%%%%%%%%%%%
persistent Alpha_pos Beta_pos Delta_pos a A1 A2 A3 C1 C2 C3 X1 X2 X3 D_alpha D_beta D_delta Flag4ub Flag4lb Nexi NexiOld f Stab mean_power Reg c Cycle  MaxCycle Pbest  pass fin test  PN dc dc_Fitness Alpha_score Beta_score Delta_score w c1 c2 r1 r2 Xmin Xmax D Partical_Best_Fitness Global_best_position Global_best_fitness;

dataType = 'double';


if ACT==0 || fin==1 % Re-initialization for another test
    Nexi=0;
    NexiOld=0;
    f=0;
    Reg=1;
    Cycle=1;
    Pbest=zeros(1,5);
    pass=0;
    test=0;
    mean_power=-inf;
    c=0;
    
end 

if ACT==0 % First execution only
    clc
    save=1;
    robustness=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Algorithm initialization %%%%%%%%%%%%%%%%%%%%%%%%
if ACT==0  || fin==1


Stab=1;

MaxCycle=30; % Maximum number of cycle
PN=6; % Particules number

dc=zeros(1,PN);
dc_Fitness=zeros(1,PN);

Alpha_score=0;
Beta_score=0;
Delta_score=0;
Alpha_pos=0;
Beta_pos=0;
Delta_pos=0;

% w=0.4; % Inertia weight
% c1=1.6; % Personal weight
% c2=1.2; % Social weight

r1=rand(); % Random vectors [0,1]
r2=rand(); %

Xmin=0.05;  % search range of solutions minimum duty cycle and Maximum duty cycle 
Xmax=0.95; %

dc=Xmin + rand(1,PN)*(Xmax-Xmin); % Initialization of particles positions 

fin=0; % Disable re-initialization
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%****************************** Main program ******************************%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% control period regulation %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nexi=Nexi+1;
if NexiOld==0
    f=1;
else
    f=0;
end
if Nexi== (NexiOld+Stab+f)    
    NexiOld=Nexi;
    Reg=Reg+1;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Initial Population (first cycle) %%%%%%%%%%%%%%%%%%%%%%
if Reg<=PN && Cycle==1
    D=dc(1,Reg);
end
if Reg>1 && Reg<=PN+1 && Nexi== NexiOld && Cycle==1 
    dc_Fitness(1,Reg-1)=P;
end
if Reg==PN+1 && Nexi== NexiOld && Cycle==1
  for i=1:PN          
        % Update Alpha, Beta, and Delta
        if dc_Fitness(1,i)>Alpha_score 
            Alpha_score=dc_Fitness(1,i); % Update alpha
            Alpha_pos=dc(1,i);
        end
        
        if dc_Fitness(1,i)<Alpha_score && dc_Fitness(1,i)>Beta_score 
            Beta_score=dc_Fitness(1,i); % Update beta
            Beta_pos=dc(1,i);
        end
        
        if dc_Fitness(1,i)<Alpha_score && dc_Fitness(1,i)<Beta_score && dc_Fitness(1,i)>Delta_score 
            Delta_score=dc_Fitness(1,i); % Update delta
            Delta_pos=dc(1,i);
        end
 end

 % Update the Position of search agents including omegas
         a=2-Cycle*((2)/MaxCycle); % a decreases linearly fron 2 to 0
 for i=1:PN  
           
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % Equation (3.3)
            C1=2*r2; % Equation (3.4)
            
            D_alpha=abs(C1*Alpha_pos-dc(1,i)); % Equation (3.5)-part 1
            X1=Alpha_pos-A1*D_alpha; % Equation (3.6)-part 1
                       
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2; % Equation (3.4)
            
            D_beta=abs(C2*Beta_pos-dc(1,i)); % Equation (3.5)-part 2
            X2=Beta_pos-A2*D_beta; % Equation (3.6)-part 2       
            
            r1=rand();
            r2=rand(); 
            
            A3=2*a*r1-a; % Equation (3.3)
            C3=2*r2; % Equation (3.4)
            
            D_delta=abs(C3*Delta_pos-dc(1,i)); % Equation (3.5)-part 3
            X3=Delta_pos-A3*D_delta; % Equation (3.5)-part 3             
            
            dc(1,i)=(X1+X2+X3)/3;% Equation (3.7)
            
                  Flag4ub=dc(1,i)>Xmax; 
                  Flag4lb=dc(1,i)<Xmin;
                  dc(1,i)=(dc(1,i).*(~(Flag4ub+Flag4lb)))+Xmax.*Flag4ub+Xmin.*Flag4lb;               
    
 end
            
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Other cycles %%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Cycle<=MaxCycle
    
if Reg>=PN+1 && Reg<=2*PN
    D=dc(1,Reg-PN);
end
if Reg>PN+1 && Reg<=(2*PN)+1 && Nexi== NexiOld  
    dc_Fitness(1,Reg-PN-1)=P;
end

if Reg==(2*PN)+1 && Nexi== NexiOld
for i=1:PN          
        % Update Alpha, Beta, and Delta
        if dc_Fitness(1,i)>Alpha_score 
            Alpha_score=dc_Fitness(1,i); % Update alpha
            Alpha_pos=dc(1,i);
        end
        
        if dc_Fitness(1,i)<Alpha_score && dc_Fitness(1,i)>Beta_score 
            Beta_score=dc_Fitness(1,i); % Update beta
            Beta_pos=dc(1,i);
        end
        
        if dc_Fitness(1,i)<Alpha_score && dc_Fitness(1,i)<Beta_score && dc_Fitness(1,i)>Delta_score 
            Delta_score=dc_Fitness(1,i); % Update delta
            Delta_pos=dc(1,i);
        end
end
% end
% 
% if Reg==(2*PN)+1 && Nexi== NexiOld
    a=2-Cycle*((2)/MaxCycle);
for i=1:PN  
          % a decreases linearly fron 2 to 0
           
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % Equation (3.3)
            C1=2*r2; % Equation (3.4)
            
            D_alpha=abs(C1*Alpha_pos-dc(1,i)); % Equation (3.5)-part 1
            X1=Alpha_pos-A1*D_alpha; % Equation (3.6)-part 1
                       
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2; % Equation (3.4)
            
            D_beta=abs(C2*Beta_pos-dc(1,i)); % Equation (3.5)-part 2
            X2=Beta_pos-A2*D_beta; % Equation (3.6)-part 2       
            
            r1=rand();
            r2=rand(); 
            
            A3=2*a*r1-a; % Equation (3.3)
            C3=2*r2; % Equation (3.4)
            
            D_delta=abs(C3*Delta_pos-dc(1,i)); % Equation (3.5)-part 3
            X3=Delta_pos-A3*D_delta; % Equation (3.5)-part 3             
            
            dc(1,i)=(X1+X2+X3)/3;% Equation (3.7)
            
                  Flag4ub=dc(1,i)>Xmax; 
                  Flag4lb=dc(1,i)<Xmin;
                  dc(1,i)=(dc(1,i).*(~(Flag4ub+Flag4lb)))+Xmax.*Flag4ub+Xmin.*Flag4lb;               
    
 end

Reg=PN;
D=Alpha_pos;

Cycle=Cycle+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% stop criterion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% STD of last 7 achieved maximum power (7 last cycles)

Pbest(1,1:4)=Pbest(1,2:5);
Pbest(1,5)=Alpha_score;
STD=std(Pbest);
mean_power=mean(Pbest);

if (Cycle>=MaxCycle || STD<=0.5) 
    pass=0;
    Cycle=MaxCycle+1;
    D=Alpha_pos;
    STD=100;
    c=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
end

%%%%%%%%%%%%%%%% Re-initialization for another test %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


test=(P-mean_power)/mean_power;

if (test>0.1 || test<-0.1) && c==1 && pass>=5
    fin=1;
end

pass=pass+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





out=D; % Output duty cycle



end