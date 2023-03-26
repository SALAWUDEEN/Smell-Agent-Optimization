%___________________________________________________________________________%
%  Smell Agent Optimization (SAO) source codes Verion v1.0                   %
%                                                                           %
%  Developed in MATLAB R2017a                                               %
%                                                                           %
%  Author and programmer: Salawudeen Ahmed Tijani                           %
%                                                                           %
%         e-Mail: tasalawudeen@abu.edu.ng                                   %
%                                                                           %
%  Credits: Mr. Salawudeen Ahmed Tijani                                     %
%           Prof. Muhammad Bashir Mua'azu                                   %
%           Dr. Yusuf Abubakar Sha'aban                                     %
%           Dr. Adedokun Emmanuel Adewale                                   %                                                                           %
%__________________________________________________________________________
%_%
function SAO_run
clc;close all;
format long

N=50;%input('Provide the Number of Smell Molecules:= ');%INITIAL POPULATION
clc
T=3;%Temperature of gas molecules.

K=1.38064852*10^(-23);% This is the boltzman constant
run=1000;
SN=2.5;
Fname='F10'; %Select function
[lb,ub,dim,myCost] = Select_Function(Fname);
m=2.4;%This is the mass of the molecules.
D=dim;
% Initial population (position) of the gas molecules as follows
BestCost=inf;
tic
for k=1:run    
    
    for i=1:N
    for j=1:D               
        molecules(i,j)=lb(j)+rand()*(ub(j)-lb(j));%Defind the initial positions of the smell molecules.    
    end
    end
% molecules=rand(N,D);
v=molecules*0.1;
molecules=molecules+v;%Initial populaation of SAO
for i=1:N
    y(i)=myCost(molecules(i,:));%Evaluate the fitness of the initial smell molecules
    if y(i)<BestCost
        BestCost=y(1);
    end
end
[ymin,index]=min(y);%Obtain the fitness of the best molecule
x_agent=molecules(index,:);%Determin the agent
olf=3.5;%Olfaction capacity of the agent.

    for i=1:N
        for j=1:D
%             Update the molecular Velocity
            v(i,j)=v(i,j)+rand*sqrt(3*K*T/m);
        end
    end 
%     perform sniffing
    for i=1:N
        for j=1:D
            molecules(i,j)=molecules(i,j)+v(i,j);
        end
    end   
    for i=1:N
        ys(i)=myCost(molecules(i,:));
        if ys(i)<BestCost
            BestCost=ys(i);
        end
    end
    [ysmin,sindex]=min(ys);
    xs_agent=molecules(sindex,:);
    [ysmax,sidx]=max(y);
    x_worst=molecules(sidx,:);%Determin the position of worst smell molecule
    if ysmin<ymin
        xs_agent=molecules(sindex,:);
        ymin=ysmin;
    end
    %%
    %Evaluate the Trailing mode
    for i=1:N
        for j=i:D
            molecules(i,j)=molecules(i,j)+rand*olf*(x_agent(1,j)-abs(molecules(i,j)))...
                -rand*olf*(x_worst(1,j)-abs(molecules(i,j)));
        end
    end
    %Make sure no smell molecules excape the boundary 
    for i=1:N
        for j=1:D
            if molecules(i,j)<lb(j)
                molecules(i,j)=lb(j);          
            elseif molecules(i,j)>ub(j)
                molecules(i,j)=ub(j);
            end
        end
    end
    %Evaluate the fitness of the Trailing mode
for i=1:N
      yt(i)=myCost(molecules(i,:));
      if yt(i)<BestCost
          BestCost=yt(i);
      end
end
    %%
%     Compare the fitness of the trailing mode and the sniffing mode
%   and implement the random mode    
    for i=1:N
        for j=1:D
            if yt(i) < ys(i)
                Best_Molecule(i,j)=BestCost;
            elseif yt(i) > ys(i)               
                molecules(i,j)=molecules(i,j)+rand()*SN;
                molecules(i,j)=molecules(i,j)+(v(i,j)+rand*sqrt(3*K*T/m));
                molecules(i,j)=molecules(i,j)+rand*olf*(x_agent(1,j)-abs(molecules(i,j)))...
                -rand*olf*(x_worst(1,j)-abs(molecules(i,j)));
            end
        end
    end
    for i=1:N
        ybest(i)=myCost(molecules(i,:));
        if ybest(i)<BestCost
            BestCost=ybest(i);           
        end
    end
    [SmellObject,Position]=min(BestCost);   
%     if iteration==1

Object(k)=sort(SmellObject,'descend');
disp(['Iteration ' num2str(k) ': Smell Object = ' num2str(Object(k))]);

end

toc
disp('The Smell Object is Obtained as: ');
Best_Object=sort(Object,'descend')';
disp('The Object Evapourating the Smell is')
Smell_Object=min(Best_Object)

figure(1)
plot(Best_Object,'k','LineWidth',2);
title('Optimization process','fontsize',12)
xlabel('Iteration Number','fontsize',12);ylabel('Smell','fontsize',12);
figure(2)
semilogy(Best_Object,'b','LineWidth',2)
title('Optimization process','fontsize',12)
xlabel('Iteration Number','fontsize',12);ylabel('Smell','fontsize',12);

