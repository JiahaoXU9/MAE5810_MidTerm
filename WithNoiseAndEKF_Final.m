%%%%clear the environment
clc;clear;
clf;
close all;

% One case to consider: can add more case
scenario_type='baseline';

[Uacc,Uomega]=get_controlinputs(scenario_type);

%state vector: x,y 2D position, velocity, heading
%simulate no noise system
[Xnonoise,n,t,dt,nt]=simulate_2Dcar(Uacc,Uomega);

%plots
plot_birdseyeview(Xnonoise,[],[],'Truth: Birds Eye View');
ii_plot=[3 4];
plot_estimator(t,Xnonoise,[],[],ii_plot,'Truth: Velocity/Heading States');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create accel/rate gyro measurements with bias = 0
bias_acc=0.0;
bias_rg=0.0;
nw=2;
Qw=diag([0.1^2 0.1^2]);
w=sqrtm(Qw)*randn(nw,nt);
Zacc=Uacc+bias_acc+w(1,:);
Zrg=Uomega+bias_rg+w(2,:);

%Create 2D GPS like measurements 
nz=2;
R=eye(nz)*0.1^2;
v=sqrtm(R)*randn(nz,nt);
Z=[Xnonoise(1:2,:)] + v;
Hgps=[eye(2) zeros(2,2)]; %output matrix is linear


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extended Kalman Filter (EKF)
x0=[0;0;0;0]; %no bias case
P0=diag([2^2 2^2 1^2 0.1^2]); %no bias case
n=length(x0);
Q=Qw; 
H=Hgps;
xhatp=x0;Pp(1:n,1:n,1)=P0;
xhatu=x0;Pu(1:n,1:n,1)=P0;
for k=1:(nt-1),    
    %predict state
    xhatp(:,k+1)=predict_state_carpose(xhatu(:,k),[Zacc(k);Zrg(k)],dt);
    %predict covariance    
    [F,G]=getFG_carpose(xhatu(:,k),dt);    
    Pp(1:n,1:n,k+1) = F*Pu(1:n,1:n,k)*F' + G*Q*G';
    %Kalman Gain
    K = Pp(1:n,1:n,k+1)*H'*inv(H*Pp(1:n,1:n,k+1)*H' + R);
    %Update
    xhatu(:,k+1) = xhatp(:,k+1) + K *(Z(:,k+1) - H*xhatp(:,k+1));
    Pu(1:n,1:n,k+1) = (eye(n)-K*H)*Pp(1:n,1:n,k+1)*(eye(n)-K*H)' + K*R*K';
end
%
ii_plot=[1 2];
plot_estimator_error(t,Xnonoise,xhatu,Pu,ii_plot,'EKF: North/East States');
ii_plot=[3 4];
plot_estimator_error(t,Xnonoise,xhatu,Pu,ii_plot,'EKF: Velocity/Heading States');
plot_birdseyeview(Xnonoise,xhatu,Pu,'EKF: Birds Eye');

n=nt;
dt=0.1;
xtrue=xhatu(1,:);
ytrue=xhatu(2,:);
xvtrue=xhatu(3,:).*cos(xhatu(4,:));
yvtrue=xhatu(3,:).*sin(xhatu(4,:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%SENSOR ONE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%the position of Oa in world frame
Oaw=[2;5;8];

positionTargetRealAWorld = [xtrue' ytrue' zeros(nt,1)];
positionTargetRealAInertial = positionTargetRealAWorld-repmat(Oaw',nt,1);

xoatrue = xtrue;
yoatrue = ytrue;
for i=1:n
    if xoatrue(i)>7
        xoatrue(i) = xoatrue(i-1);
        yoatrue(i) = yoatrue(i-1);
    end
end


        

%%%%the position of 4 corners of the virtual image plane in the body frame
q1b=[0.018;0.0135;0.05];
q2b=[0.018;-0.0135;0.05];
q3b=[-0.018;-0.0135;0.05];
q4b=[-0.018;0.0135;0.05];
centerb =[0;0;0.05];

%%%%calcuate thr ideal fai and kai
apositionTargetWorld = [xoatrue' yoatrue' zeros(nt,1)];
apositionTargetInertial = apositionTargetWorld-repmat(Oaw',nt,1);%now it represent a vector from the pinhole to the target position
apositionTargetUnitInertial = apositionTargetInertial./sqrt(apositionTargetInertial(:,1).^2+apositionTargetInertial(:,2).^2+apositionTargetInertial(:,3).^2);
afaiTarget = zeros(nt,1);
afaiTarget = acosd(apositionTargetUnitInertial(:,3));
akaiTarget = zeros(nt,1);
akaiTarget=asind(apositionTargetUnitInertial(:,1)./sind(afaiTarget));

%%%%we use LQR so we need to caculate K using a discrete system
Ad=[1 0 dt 0;
    0 1 0 dt;
    0 0 1 0;
    0 0 0 1;];
Bd = [0 0;
    0 0;
    1.7453   0;
   0 1.7453  ;];
Q = [0.00001 0 0 0;
    0 0.000001 0 0;
    0 0 0 0;
    0 0 0 0;];
R = [100000 0;
    0 100000;];
aK =lqrd(Ad,Bd,Q,R,dt);

%%%%do some space initializing
afai = zeros(1,n);
akai = zeros(1,n);
afaidot = zeros(1,n);
akaidot = zeros(1,n);
au1 = zeros(1,n);
au2 = zeros(1,n);
aq1W=zeros(3,n);
aq2W=zeros(3,n);
aq3W=zeros(3,n);
aq4W=zeros(3,n);
acenterW = zeros(3,n);
afai(1) = afaiTarget(1);
akai(1) = akaiTarget(1);
afaidot(1) = 0;
akaidot(1) = 0;
au1(1) = 0;
au2(1) = 0;
step = 1;

%%%%calculate the FOV of the first time-step
aq1new = funkai(akai(1))*funfai(afai(1))*q1b;
aq2new = funkai(akai(1))*funfai(afai(1))*q2b;
aq3new = funkai(akai(1))*funfai(afai(1))*q3b;
aq4new = funkai(akai(1))*funfai(afai(1))*q4b;
acenternew = funkai(akai(1))*funfai(afai(1))*centerb;

aq1original = [aq1new(1)+Oaw(1);aq1new(2)+Oaw(2);aq1new(3)+Oaw(3)];
aq2original = [aq2new(1)+Oaw(1);aq2new(2)+Oaw(2);aq2new(3)+Oaw(3)];
aq3original = [aq3new(1)+Oaw(1);aq3new(2)+Oaw(2);aq3new(3)+Oaw(3)];
aq4original = [aq4new(1)+Oaw(1);aq4new(2)+Oaw(2);aq4new(3)+Oaw(3)];
acenteroriginal = [acenternew(1)+Oaw(1);acenternew(2)+Oaw(2);acenternew(3)+Oaw(3)];

aq1W(:,1) = [(-aq1original(3))/(Oaw(3)-aq1original(3))*(Oaw(1)-aq1original(1))+aq1original(1);(-aq1original(3))/(Oaw(3)-aq1original(3))*(Oaw(2)-aq1original(2))+aq1original(2);0];
aq2W(:,1) = [(-aq2original(3))/(Oaw(3)-aq2original(3))*(Oaw(1)-aq2original(1))+aq2original(1);(-aq2original(3))/(Oaw(3)-aq2original(3))*(Oaw(2)-aq2original(2))+aq2original(2);0];
aq3W(:,1) = [(-aq3original(3))/(Oaw(3)-aq3original(3))*(Oaw(1)-aq3original(1))+aq3original(1);(-aq3original(3))/(Oaw(3)-aq3original(3))*(Oaw(2)-aq3original(2))+aq3original(2);0];
aq4W(:,1) = [(-aq4original(3))/(Oaw(3)-aq4original(3))*(Oaw(1)-aq4original(1))+aq4original(1);(-aq4original(3))/(Oaw(3)-aq4original(3))*(Oaw(2)-aq4original(2))+aq4original(2);0];
acenterW(:,1) = [(-acenteroriginal(3))/(Oaw(3)-acenteroriginal(3))*(Oaw(1)-acenteroriginal(1))+acenteroriginal(1);(-acenteroriginal(3))/(Oaw(3)-acenteroriginal(3))*(Oaw(2)-acenteroriginal(2))+acenteroriginal(2);0];

%%%%doing the whole simulation
for i =0.1:0.1:n/10-0.1
    akai(step+1) =akai(step)+dt*akaidot(step);
    afai(step+1)=afai(step)+dt*afaidot(step);
    akaidot(step+1)=akaidot(step)+au1(step)*dt;
    afaidot(step+1)=afaidot(step)+au2(step)*dt;
    akaiideal = akaiTarget(step+1);
    afaiideal = afaiTarget(step+1);
    adeltakai = akai(step+1)-akaiideal;
    adeltafai = afai(step+1)-afaiideal;
    aresult = -aK*[adeltakai;adeltafai;akaidot(step+1);afaidot(step+1)];
    au1(step+1) = aresult(1);
    au2(step+1) = aresult(2);
    
    %%%%caculate the FOV
    aq1new = funkai(akai(step+1))*funfai(afai(step+1))*q1b;
    aq2new = funkai(akai(step+1))*funfai(afai(step+1))*q2b;
    aq3new = funkai(akai(step+1))*funfai(afai(step+1))*q3b;
    aq4new = funkai(akai(step+1))*funfai(afai(step+1))*q4b;
    acenternew = funkai(akai(step+1))*funfai(afai(step+1))*centerb;

    aq1original = [aq1new(1)+Oaw(1);aq1new(2)+Oaw(2);aq1new(3)+Oaw(3)];
    aq2original = [aq2new(1)+Oaw(1);aq2new(2)+Oaw(2);aq2new(3)+Oaw(3)];
    aq3original = [aq3new(1)+Oaw(1);aq3new(2)+Oaw(2);aq3new(3)+Oaw(3)];
    aq4original = [aq4new(1)+Oaw(1);aq4new(2)+Oaw(2);aq4new(3)+Oaw(3)];
    acenteroriginal = [acenternew(1)+Oaw(1);acenternew(2)+Oaw(2);acenternew(3)+Oaw(3)];

    aq1W(:,step+1) = [(-aq1original(3))/(Oaw(3)-aq1original(3))*(Oaw(1)-aq1original(1))+aq1original(1);(-aq1original(3))/(Oaw(3)-aq1original(3))*(Oaw(2)-aq1original(2))+aq1original(2);0];
    aq2W(:,step+1) = [(-aq2original(3))/(Oaw(3)-aq2original(3))*(Oaw(1)-aq2original(1))+aq2original(1);(-aq2original(3))/(Oaw(3)-aq2original(3))*(Oaw(2)-aq2original(2))+aq2original(2);0];
    aq3W(:,step+1) = [(-aq3original(3))/(Oaw(3)-aq3original(3))*(Oaw(1)-aq3original(1))+aq3original(1);(-aq3original(3))/(Oaw(3)-aq3original(3))*(Oaw(2)-aq3original(2))+aq3original(2);0];
    aq4W(:,step+1) = [(-aq4original(3))/(Oaw(3)-aq4original(3))*(Oaw(1)-aq4original(1))+aq4original(1);(-aq4original(3))/(Oaw(3)-aq4original(3))*(Oaw(2)-aq4original(2))+aq4original(2);0];
    acenterW(:,step+1) = [(-acenteroriginal(3))/(Oaw(3)-acenteroriginal(3))*(Oaw(1)-acenteroriginal(1))+acenteroriginal(1);(-acenteroriginal(3))/(Oaw(3)-acenteroriginal(3))*(Oaw(2)-acenteroriginal(2))+acenteroriginal(2);0];
    
    step = step+1;

end


%%%doing some ploting
plotme(akai,akaiTarget,'tilt','tilttarget',n)
ylabel("Tilt")
title("Camera 1 Tilt")

plotme(afai,afaiTarget,'pan','pantarget',n)
ylabel("Pan")
title("Camera 1 Pan")

plotme(akaidot,afaidot,'kaidot','faidot',n)
ylabel("angular velocity")
title("Camera 1 angular velocity")

plotme(au1,au2,'u1','u2',n)
ylabel("input")
title("Camera 1 control input")



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%SENSOR TWO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%the position of Oa in world frame
Obw=[12;9;8];

positionTargetRealWorld = [xtrue' ytrue' zeros(nt,1)];
positionTargetRealInertial = positionTargetRealWorld-repmat(Obw',nt,1);

xobtrue = xtrue;
yobtrue = ytrue;
for i=n-1:-1:1
    if xobtrue(i)<7
        xobtrue(i) = xobtrue(i+1);
        yobtrue(i) = yobtrue(i+1);
    end
end


        

%%%%the position of 4 corners of the virtual image plane in the body frame
q1b=[0.018;0.0135;0.05];
q2b=[0.018;-0.0135;0.05];
q3b=[-0.018;-0.0135;0.05];
q4b=[-0.018;0.0135;0.05];
centerb =[0;0;0.05];

%%%%calcuate thr ideal fai and kai
positionTargetWorld = [xobtrue' yobtrue' zeros(nt,1)];
positionTargetInertial = positionTargetWorld-repmat(Obw',nt,1);%now it represent a vector from the pinhole to the target position
positionTargetUnitInertial = positionTargetInertial./sqrt(positionTargetInertial(:,1).^2+positionTargetInertial(:,2).^2+positionTargetInertial(:,3).^2);
faiTarget = zeros(nt,1);
faiTarget = acosd(positionTargetUnitInertial(:,3));
kaiTarget = zeros(nt,1);
kaiTarget=asind(positionTargetUnitInertial(:,1)./sind(faiTarget));

%%%%we use LQR so we need to caculate K using a discrete system
Ad=[1 0 dt 0;
    0 1 0 dt;
    0 0 1 0;
    0 0 0 1;];
Bd = [0 0;
    0 0;
    1.7453   0;
   0 1.7453  ;];
Q = [0.00001 0 0 0;
    0 0.000001 0 0;
    0 0 0 0;
    0 0 0 0;];
R = [100000 0;
    0 100000;];
K =lqrd(Ad,Bd,Q,R,dt);

%%%%do some space initializing
fai = zeros(1,n);
kai = zeros(1,n);
faidot = zeros(1,n);
kaidot = zeros(1,n);
u1 = zeros(1,n);
u2 = zeros(1,n);
q1W=zeros(3,n);
q2W=zeros(3,n);
q3W=zeros(3,n);
q4W=zeros(3,n);
centerW = zeros(3,n);
fai(1) = faiTarget(1);
kai(1) = kaiTarget(1);
faidot(1) = 0;
kaidot(1) = 0;
u1(1) = 0;
u2(1) = 0;
step = 1;

%%%%calculate the FOV of the first time-step
q1new = funkai(kai(1))*funfai(fai(1))*q1b;
q2new = funkai(kai(1))*funfai(fai(1))*q2b;
q3new = funkai(kai(1))*funfai(fai(1))*q3b;
q4new = funkai(kai(1))*funfai(fai(1))*q4b;
centernew = funkai(kai(1))*funfai(fai(1))*centerb;

q1original = [q1new(1)+Obw(1);q1new(2)+Obw(2);q1new(3)+Obw(3)];
q2original = [q2new(1)+Obw(1);q2new(2)+Obw(2);q2new(3)+Obw(3)];
q3original = [q3new(1)+Obw(1);q3new(2)+Obw(2);q3new(3)+Obw(3)];
q4original = [q4new(1)+Obw(1);q4new(2)+Obw(2);q4new(3)+Obw(3)];
centeroriginal = [centernew(1)+Obw(1);centernew(2)+Obw(2);centernew(3)+Obw(3)];

q1W(:,1) = [(-q1original(3))/(Obw(3)-q1original(3))*(Obw(1)-q1original(1))+q1original(1);(-q1original(3))/(Obw(3)-q1original(3))*(Obw(2)-q1original(2))+q1original(2);0];
q2W(:,1) = [(-q2original(3))/(Obw(3)-q2original(3))*(Obw(1)-q2original(1))+q2original(1);(-q2original(3))/(Obw(3)-q2original(3))*(Obw(2)-q2original(2))+q2original(2);0];
q3W(:,1) = [(-q3original(3))/(Obw(3)-q3original(3))*(Obw(1)-q3original(1))+q3original(1);(-q3original(3))/(Obw(3)-q3original(3))*(Obw(2)-q3original(2))+q3original(2);0];
q4W(:,1) = [(-q4original(3))/(Obw(3)-q4original(3))*(Obw(1)-q4original(1))+q4original(1);(-q4original(3))/(Obw(3)-q4original(3))*(Obw(2)-q4original(2))+q4original(2);0];
centerW(:,1) = [(-centeroriginal(3))/(Obw(3)-centeroriginal(3))*(Obw(1)-centeroriginal(1))+centeroriginal(1);(-centeroriginal(3))/(Obw(3)-centeroriginal(3))*(Obw(2)-centeroriginal(2))+centeroriginal(2);0];

%%%%doing the whole simulation
for i =0.1:0.1:n/10-0.1
    kai(step+1) =kai(step)+dt*kaidot(step);
    fai(step+1)=fai(step)+dt*faidot(step);
    kaidot(step+1)=kaidot(step)+u1(step)*dt;
    faidot(step+1)=faidot(step)+u2(step)*dt;
    kaiideal = kaiTarget(step+1);
    faiideal = faiTarget(step+1);
    deltakai = kai(step+1)-kaiideal;
    deltafai = fai(step+1)-faiideal;
    result = -K*[deltakai;deltafai;kaidot(step+1);faidot(step+1)];
    u1(step+1) = result(1);
    u2(step+1) = result(2);
    
    %%%%caculate the FOV
    q1new = funkai(kai(step+1))*funfai(fai(step+1))*q1b;
    q2new = funkai(kai(step+1))*funfai(fai(step+1))*q2b;
    q3new = funkai(kai(step+1))*funfai(fai(step+1))*q3b;
    q4new = funkai(kai(step+1))*funfai(fai(step+1))*q4b;
    centernew = funkai(kai(step+1))*funfai(fai(step+1))*centerb;

    q1original = [q1new(1)+Obw(1);q1new(2)+Obw(2);q1new(3)+Obw(3)];
    q2original = [q2new(1)+Obw(1);q2new(2)+Obw(2);q2new(3)+Obw(3)];
    q3original = [q3new(1)+Obw(1);q3new(2)+Obw(2);q3new(3)+Obw(3)];
    q4original = [q4new(1)+Obw(1);q4new(2)+Obw(2);q4new(3)+Obw(3)];
    centeroriginal = [centernew(1)+Obw(1);centernew(2)+Obw(2);centernew(3)+Obw(3)];

    q1W(:,step+1) = [(-q1original(3))/(Obw(3)-q1original(3))*(Obw(1)-q1original(1))+q1original(1);(-q1original(3))/(Obw(3)-q1original(3))*(Obw(2)-q1original(2))+q1original(2);0];
    q2W(:,step+1) = [(-q2original(3))/(Obw(3)-q2original(3))*(Obw(1)-q2original(1))+q2original(1);(-q2original(3))/(Obw(3)-q2original(3))*(Obw(2)-q2original(2))+q2original(2);0];
    q3W(:,step+1) = [(-q3original(3))/(Obw(3)-q3original(3))*(Obw(1)-q3original(1))+q3original(1);(-q3original(3))/(Obw(3)-q3original(3))*(Obw(2)-q3original(2))+q3original(2);0];
    q4W(:,step+1) = [(-q4original(3))/(Obw(3)-q4original(3))*(Obw(1)-q4original(1))+q4original(1);(-q4original(3))/(Obw(3)-q4original(3))*(Obw(2)-q4original(2))+q4original(2);0];
    centerW(:,step+1) = [(-centeroriginal(3))/(Obw(3)-centeroriginal(3))*(Obw(1)-centeroriginal(1))+centeroriginal(1);(-centeroriginal(3))/(Obw(3)-centeroriginal(3))*(Obw(2)-centeroriginal(2))+centeroriginal(2);0];
    
    step = step+1;

end


%%%%doing some ploting
plotme(kai,kaiTarget,'tilt','tilttarget',n)
ylabel("Tilt")
title("Camera 2 Tilt")

plotme(fai,faiTarget,'pan','pantarget',n)
ylabel("Pan")
title("Camera 2 Pan")

plotme(kaidot,faidot,'kaidot','faidot',n)
ylabel("angular velocity")
title("Camera 2 angular velocity")

plotme(u1,u2,'u1','u2',n)
ylabel("input")
title("Camera 2 control input")



% %%%%%%%%NOW DO THE ANIMATION
%Total 700 frame
clc
figure()
M = moviein(n);

%Plot everyframe

for i=1:n
    clf
    hold on
    grid on
        title("Tracking using EKF in the case with Noise and No Bias")
    view([1 1 1])
    boolean = isinfov(positionTargetRealAInertial(i,:),akai(i),afai(i));
    if xoatrue(i)>=xoatrue(700)
         text(30,-10,10,'Camera 1 Stop Tracking')
         text(30,-10,9,'Camera 2 is NOW Tracking')
    else
         text(30,-10,10,'Camera 1 is NOW Tracking')
    end
    text(1,5,8,'Camera 1')
    text(12,3,8,'Camera 2')
     
        
        
    if boolean == 0        
%         fprintf("Target not in Field-of-View!!")
        patch([aq1W(1,i),aq2W(1,i),aq3W(1,i),aq4W(1,i)],[aq1W(2,i),aq2W(2,i),aq3W(2,i),aq4W(2,i)],'red','FaceAlpha',.3)
        plot3([aq1W(1,i),aq2W(1,i)],[aq1W(2,i),aq2W(2,i)],[aq1W(3,i),aq2W(3,i)],'red','LineWidth',1)
        plot3([aq2W(1,i),aq3W(1,i)],[aq2W(2,i),aq3W(2,i)],[aq2W(3,i),aq3W(3,i)],'red','LineWidth',1)
        plot3([aq3W(1,i),aq4W(1,i)],[aq3W(2,i),aq4W(2,i)],[aq3W(3,i),aq4W(3,i)],'red','LineWidth',1)
        plot3([aq4W(1,i),aq1W(1,i)],[aq4W(2,i),aq1W(2,i)],[aq4W(3,i),aq1W(3,i)],'red','LineWidth',1)
    else
        patch([aq1W(1,i),aq2W(1,i),aq3W(1,i),aq4W(1,i)],[aq1W(2,i),aq2W(2,i),aq3W(2,i),aq4W(2,i)],'green','FaceAlpha',.3)
        plot3([aq1W(1,i),aq2W(1,i)],[aq1W(2,i),aq2W(2,i)],[aq1W(3,i),aq2W(3,i)],'g','LineWidth',1)
        plot3([aq2W(1,i),aq3W(1,i)],[aq2W(2,i),aq3W(2,i)],[aq2W(3,i),aq3W(3,i)],'g','LineWidth',1)
        plot3([aq3W(1,i),aq4W(1,i)],[aq3W(2,i),aq4W(2,i)],[aq3W(3,i),aq4W(3,i)],'g','LineWidth',1)
        plot3([aq4W(1,i),aq1W(1,i)],[aq4W(2,i),aq1W(2,i)],[aq4W(3,i),aq1W(3,i)],'g','LineWidth',1)
    end
    
    plot3([Xnonoise(1,i)],[Xnonoise(2,i)],[0],'r.','markersize',25)
    plot3([Oaw(1)],[Oaw(2)],[Oaw(3)],'b.','markersize',25)
    plot3([aq1W(1,i)],[aq1W(2,i)],[aq1W(3,i)],'k.','markersize',10)
    plot3([aq2W(1,i)],[aq2W(2,i)],[aq2W(3,i)],'k.','markersize',10)
    plot3([aq3W(1,i)],[aq3W(2,i)],[aq3W(3,i)],'k.','markersize',10)
    plot3([aq4W(1,i)],[aq4W(2,i)],[aq4W(3,i)],'k.','markersize',10)
    
    plot3([Oaw(1),aq1W(1,i)],[Oaw(2),aq1W(2,i)],[Oaw(3),aq1W(3,i)],'k-.','LineWidth',1)
    plot3([Oaw(1),aq2W(1,i)],[Oaw(2),aq2W(2,i)],[Oaw(3),aq2W(3,i)],'k-.','LineWidth',1)
    plot3([Oaw(1),aq3W(1,i)],[Oaw(2),aq3W(2,i)],[Oaw(3),aq3W(3,i)],'k-.','LineWidth',1)
    plot3([Oaw(1),aq4W(1,i)],[Oaw(2),aq4W(2,i)],[Oaw(3),aq4W(3,i)],'k-.','LineWidth',1)
    plot3([Oaw(1),acenterW(1,i)],[Oaw(2),acenterW(2,i)],[Oaw(3),acenterW(3,i)],'b','LineWidth',1)
    
    
    
    
    
     boolean = isinfov(positionTargetRealInertial(i,:),kai(i),fai(i));
    if boolean == 0        
%         fprintf("Target not in Field-of-View!!")
        patch([q1W(1,i),q2W(1,i),q3W(1,i),q4W(1,i)],[q1W(2,i),q2W(2,i),q3W(2,i),q4W(2,i)],'red','FaceAlpha',.3)
        plot3([q1W(1,i),q2W(1,i)],[q1W(2,i),q2W(2,i)],[q1W(3,i),q2W(3,i)],'red','LineWidth',1)
        plot3([q2W(1,i),q3W(1,i)],[q2W(2,i),q3W(2,i)],[q2W(3,i),q3W(3,i)],'red','LineWidth',1)
        plot3([q3W(1,i),q4W(1,i)],[q3W(2,i),q4W(2,i)],[q3W(3,i),q4W(3,i)],'red','LineWidth',1)
        plot3([q4W(1,i),q1W(1,i)],[q4W(2,i),q1W(2,i)],[q4W(3,i),q1W(3,i)],'red','LineWidth',1)
    else
        patch([q1W(1,i),q2W(1,i),q3W(1,i),q4W(1,i)],[q1W(2,i),q2W(2,i),q3W(2,i),q4W(2,i)],'green','FaceAlpha',.3)
        plot3([q1W(1,i),q2W(1,i)],[q1W(2,i),q2W(2,i)],[q1W(3,i),q2W(3,i)],'g','LineWidth',1)
        plot3([q2W(1,i),q3W(1,i)],[q2W(2,i),q3W(2,i)],[q2W(3,i),q3W(3,i)],'g','LineWidth',1)
        plot3([q3W(1,i),q4W(1,i)],[q3W(2,i),q4W(2,i)],[q3W(3,i),q4W(3,i)],'g','LineWidth',1)
        plot3([q4W(1,i),q1W(1,i)],[q4W(2,i),q1W(2,i)],[q4W(3,i),q1W(3,i)],'g','LineWidth',1)
    end
    
%     plot3([Xnonoise(1,i)],[Xnonoise(2,i)],[0],'r.','markersize',25)
    plot3([Obw(1)],[Obw(2)],[Obw(3)],'b.','markersize',25)
    plot3([q1W(1,i)],[q1W(2,i)],[q1W(3,i)],'k.','markersize',10)
    plot3([q2W(1,i)],[q2W(2,i)],[q2W(3,i)],'k.','markersize',10)
    plot3([q3W(1,i)],[q3W(2,i)],[q3W(3,i)],'k.','markersize',10)
    plot3([q4W(1,i)],[q4W(2,i)],[q4W(3,i)],'k.','markersize',10)
    
    plot3([6.97 6.97],[-15,15],[0,0],'k','LineWidth',2)
    plot3([Obw(1),q1W(1,i)],[Obw(2),q1W(2,i)],[Obw(3),q1W(3,i)],'k-.','LineWidth',1)
    plot3([Obw(1),q2W(1,i)],[Obw(2),q2W(2,i)],[Obw(3),q2W(3,i)],'k-.','LineWidth',1)
    plot3([Obw(1),q3W(1,i)],[Obw(2),q3W(2,i)],[Obw(3),q3W(3,i)],'k-.','LineWidth',1)
    plot3([Obw(1),q4W(1,i)],[Obw(2),q4W(2,i)],[Obw(3),q4W(3,i)],'k-.','LineWidth',1)
    plot3([Obw(1),centerW(1,i)],[Obw(2),centerW(2,i)],[Obw(3),centerW(3,i)],'b','LineWidth',1)
    
    
    
    
    
    
    
    
    
    
    axis([-5 25 -15 15 0 Obw(3)])
    %save image to struct M
    M(i) = getframe;
end

%Do the animation and save as MP4
writerObj=VideoWriter('ProjectWithNoiseAndEKF.mp4','MPEG-4');
open(writerObj);
for i=1:n
    writeVideo(writerObj,M(i));
end
close(writerObj);


%%%%END OF ANIMATION
















































%%%%FUNCTIONS
function boolean = isinfov(positionTargetIntertial,kai,fai)
vectorbody = funfai(fai)'*funkai(kai)'*positionTargetIntertial';
ximageplane = 0.05*vectorbody(1)/vectorbody(3);
yimageplane=0.05*vectorbody(2)/vectorbody(3);
if abs(ximageplane)<=0.018 && abs(yimageplane)<=0.0135
    boolean =  1;
else
%     fprintf("x=%d,y=%d",ximageplane,yimageplane);
    boolean = 0;
end
end

function Hfai = funfai(fai)
Hfai=[1,0,0;
    0,cosd(fai),-sind(fai);
    0,sind(fai),cosd(fai);];
end

function Hkai = funkai(kai)
Hkai=[cosd(kai),-sind(kai),0;
     sind(kai),cosd(kai),0;
     0,0,1];
end

function plotme(first,second,firstname,secondname,n);
figure()
plot(0.1:0.1:n/10,first,'b','LineWidth',2,'DisplayName',firstname);
hold on;
plot(0.1:0.1:n/10,second,'r','LineWidth',2,'DisplayName',secondname);
xlabel("Time")
legend()
end

function [Uacc,Uomega]=get_controlinputs(scenario_type);
%
%scenario_type='baseline';
%scenario_type='swervy';
%
%   generates clean acceleration and rotational rate control inputs
%
%
%define inputs: acceleration and heading rate
if strcmp(scenario_type,'baseline'),
    Uacc=[ones(1,60)*0.1 zeros(1,100) -ones(1,60)*0.1 zeros(1,20) ones(1,60)*0.07 zeros(1,50) -ones(1,60)*0.07 zeros(1,50) ones(1,60)*0.1 zeros(1,100) -ones(1,60)*0.1 zeros(1,20)];
    Uomega=[zeros(1,240) ones(1,40)*pi/2/4 zeros(1,180) ones(1,40)*-pi/2/4 zeros(1,100) ones(1,40)*pi/2/4 zeros(1,60)];
else,
    return;
end
%
end

function [Xnonoise,n,t,dt,nt]=simulate_2Dcar(Uacc,Uomega);
%
%   simulates a 2D car using acceleration and rotational rate control inputs
%   state vector: x,y 2D position, velocity, heading
%
dt=0.1;
nt=length(Uacc);
t=[0:dt:dt*(nt-1)];
n=4;
Xnonoise=zeros(n,nt);
for k=1:(nt-1),
    Vk=Xnonoise(3,k);
    Tk=Xnonoise(4,k);    
    Xnonoise(:,k+1) = Xnonoise(:,k) +...
        dt*[Vk*cos(Tk);Vk*sin(Tk);Uacc(k);Uomega(k)];   
end
end

function plot_birdseyeview(x1,x2,P2,title_name);
% x1 is the true value or reference comparison
% x2,P2 is the estimator state and covariance
% 
ii_x1=[];ii_x2=[];ii_P2=[]; %for legend
figure;
if ~isempty(x1),
    plot(x1(1,:),x1(2,:),'color',[0 0.5 0]);ii_x1=1;
end
hold on;
if ~isempty(x2),
    plot(x2(1,:),x2(2,:),'b-');ii_x2=2;
end
if ~isempty(P2),
    iell=[2 10:10:700];
    for i=1:length(iell),
        ii=iell(i);
        [Xe,Ye] = calculateEllipseCov(x2([1 2],ii),P2([1 2],[1 2],ii),3);
        plot(Xe,Ye,'m-');
    end
    ii_P2=3;
end
xlabel('North (m)');ylabel('East (m)');grid;
hold off;
legend_names={'true trajectory','estimated trajectory','3\sigma bound'};
legend(legend_names{ii_x1},legend_names{ii_x2},legend_names{ii_P2},'Location','South')
%
title(title_name);
PrepFigPresentation(gcf);
end

function PrepFigPresentation(fignum);
%
% prepares a figure for presentations
%
% Fontsize: 14
% Fontweight: bold
% LineWidth: 2
% 

figure(fignum);
fig_children=get(fignum,'children'); %find all sub-plots

for i=1:length(fig_children),
    
    set(fig_children(i),'FontSize',16);
    set(fig_children(i),'FontWeight','bold');
    
    fig_children_children=get(fig_children(i),'Children');
    set(fig_children_children,'LineWidth',2);
end
end

function plot_estimator(t,x1,x2,P2,ii_plot,title_name);
% x1 is the true value or reference comparison
% x2,P2 is the estimator state and covariance
% ii_plot: 2x1 vector of which states to plot
%
axis_names={'North (m)','East (m)','Velocity (m/sec)','Heading (rad)','Accel Bias (m/sec^2)','RG Bias (rad/sec)'};
figure;subplot(122);
ii_x1=[];ii_x2=[];ii_P2=[]; %for legend
%
for i=1:length(ii_plot),
    ii=ii_plot(i);
    subplot(1,2,i);
    hold on;
    if ~isempty(x1),
        plot(t,x1(ii,:),'color',[0 0.5 0]);ii_x1=1;
    end  
    if ~isempty(x2),
        plot(t,x2(ii,:),'b-');ii_x2=2;
    end  
    if ~isempty(P2)
        plot(t,x2(ii,:)'-2*sqrt(squeeze(P2(ii,ii,:))),'b:');
        plot(t,x2(ii,:)'+2*sqrt(squeeze(P2(ii,ii,:))),'b:');ii_P2=3;
    end
    hold off
    xlabel('time (sec)');ylabel(axis_names(ii));grid;
    xlim([0 35]);set(gca,'xtick',[0:5:35]);
end
legend_names={'true state','estimate','2\sigma bound'};
legend(legend_names{ii_x1},legend_names{ii_x2},legend_names{ii_P2},'Location','South');
%
sgtitle(title_name);
PrepFigPresentation(gcf);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EKF function call
function Xkp1=predict_state_carpose(Xk,U,dt);
%
%   car pose/localization problem: NO BIAS
%   state prediction from k to k+1
%   assumes mean of process noise is zero
%
Vk=Xk(3);
Tk=Xk(4);
%
acc=U(1);
omega=U(2);
%
Xkp1 = Xk +...
    dt*[Vk*cos(Tk);
        Vk*sin(Tk);
        acc;
        omega];
%
end

function [F,G]=getFG_carpose(X,dt);
%
%   car pose/localization problem: NO BIAS
%   find the linearized system matrixes F,G
%
Vk=X(3);
Tk=X(4);
%
F=[1 0 dt*cos(Tk) -dt*Vk*sin(Tk);
   0 1 dt*sin(Tk)  dt*Vk*cos(Tk);
   0 0 1 0;
   0 0 0 1];
G=[0 0;
   0 0;
   dt 0;
   0 dt];
%
end

function plot_estimator_error(t,x1,x2,P2,ii_plot,title_name);
% x1 is the true value or reference comparison
% x2,P2 is the estimator state and covariance
% ii_plot: 2x1 vector of which states to plot
%
axis_names={'North (m)','East (m)','Velocity (m/sec)','Heading (rad)','Accel Bias (m/sec^2)','RG Bias (rad/sec)'};
figure;subplot(122);
%
for i=1:length(ii_plot),
    ii=ii_plot(i);
    subplot(1,2,i);
    err=x2(ii,:)-x1(ii,:);
    plot(t,err,'b-');
    hold on;
    if ~isempty(P2)
        plot(t,err'-2*sqrt(squeeze(P2(ii,ii,:))),'b:');
        plot(t,zeros(length(t),1),'r--');
        plot(t,err'+2*sqrt(squeeze(P2(ii,ii,:))),'b:');
    end
    hold off
    xlabel('time (sec)');ylabel(axis_names(ii));grid;
    xlim([0 35]);set(gca,'xtick',[0:5:35]);
end
legend('estimator error','2\sigma bound','zero error','Location','South');
%
sgtitle(title_name);
PrepFigPresentation(gcf);
end


function [Xe,Ye] = calculateEllipseCov(X, P, nsig, steps) 
    %# This functions returns points to draw an ellipse 
    %# 
    %#  @param X     x,y coordinates 
    %#  @param P     covariance matrix 
    %# 
 
    error(nargchk(2, 3, nargin)); 
    if nargin<3, nsig = 1; end 
    if nargin<4, steps = 36; end 
    
    [U,S,V]=svd(P);
    s1=sqrt(S(1,1));s2=sqrt(S(2,2));angle=acos(U(1,1))*180/pi;
    x=X(1);
    y=X(2);

    %scale by nsig
    s1=nsig*s1;
    s2=nsig*s2;

    beta = angle * (pi / 180); 
    sinbeta = sin(beta); 
    cosbeta = cos(beta); 
 
    alpha = linspace(0, 360, steps)' .* (pi / 180); 
    sinalpha = sin(alpha); 
    cosalpha = cos(alpha); 
 
    Xe = x + (s1 * cosalpha * cosbeta - s2 * sinalpha * sinbeta); 
    Ye = y + (s1 * cosalpha * sinbeta + s2 * sinalpha * cosbeta); 
 
end 















































