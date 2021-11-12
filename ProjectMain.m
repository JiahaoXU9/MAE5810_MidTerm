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


n=nt;
dt=0.1
xtrue=Xnonoise(1,:);
ytrue=Xnonoise(2,:);
xvtrue=Xnonoise(3,:).*cos(Xnonoise(4,:));
yvtrue=Xnonoise(3,:).*sin(Xnonoise(4,:))


%%%%the position of Oa in world frame
Oaw=[2;5;5];

%%%%the position of 4 corners of the virtual image plane in the body frame
q1b=[0.018;0.0135;0.05];
q2b=[0.018;-0.0135;0.05];
q3b=[-0.018;-0.0135;0.05];
q4b=[-0.018;0.0135;0.05];
centerb =[0;0;0.05]

%%%%calcuate thr ideal fai and kai
positionTargetWorld = [xtrue' ytrue' zeros(nt,1)];
positionTargetInertial = positionTargetWorld-repmat(Oaw',nt,1);%now it represent a vector from the pinhole to the target position
positionTargetUnitInertial = positionTargetInertial./sqrt(positionTargetInertial(:,1).^2+positionTargetInertial(:,2).^2+positionTargetInertial(:,3).^2);
faiTarget = zeros(nt,1);
faiTarget = acosd(positionTargetUnitInertial(:,3))
kaiTarget = zeros(nt,1);
kaiTarget=asind(positionTargetUnitInertial(:,1)./sind(faiTarget))




%%%%we use LQR so we need to caculate K using a discrete system
Ad=[1 0 dt 0;
    0 1 0 dt;
    0 0 1 0;
    0 0 0 1;];
Bd = [0 0;
    0 0;
    1.7453  0;
   0 1.7453 ;];
Q = [0.00000001 0 0 0;
    0 0.000000001 0 0;
    0 0 0 0;
    0 0 0 0;];
R = [100 0;
    0 100;];
K =lqrd(Ad,Bd,Q,R,dt)

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

q1original = [q1new(1)+Oaw(1);q1new(2)+Oaw(2);q1new(3)+Oaw(3)];
q2original = [q2new(1)+Oaw(1);q2new(2)+Oaw(2);q2new(3)+Oaw(3)];
q3original = [q3new(1)+Oaw(1);q3new(2)+Oaw(2);q3new(3)+Oaw(3)];
q4original = [q4new(1)+Oaw(1);q4new(2)+Oaw(2);q4new(3)+Oaw(3)];
centeroriginal = [centernew(1)+Oaw(1);centernew(2)+Oaw(2);centernew(3)+Oaw(3)];

q1W(:,1) = [(-q1original(3))/(Oaw(3)-q1original(3))*(Oaw(1)-q1original(1))+q1original(1);(-q1original(3))/(Oaw(3)-q1original(3))*(Oaw(2)-q1original(2))+q1original(2);0];
q2W(:,1) = [(-q2original(3))/(Oaw(3)-q2original(3))*(Oaw(1)-q2original(1))+q2original(1);(-q2original(3))/(Oaw(3)-q2original(3))*(Oaw(2)-q2original(2))+q2original(2);0];
q3W(:,1) = [(-q3original(3))/(Oaw(3)-q3original(3))*(Oaw(1)-q3original(1))+q3original(1);(-q3original(3))/(Oaw(3)-q3original(3))*(Oaw(2)-q3original(2))+q3original(2);0];
q4W(:,1) = [(-q4original(3))/(Oaw(3)-q4original(3))*(Oaw(1)-q4original(1))+q4original(1);(-q4original(3))/(Oaw(3)-q4original(3))*(Oaw(2)-q4original(2))+q4original(2);0];
centerW(:,1) = [(-centeroriginal(3))/(Oaw(3)-centeroriginal(3))*(Oaw(1)-centeroriginal(1))+centeroriginal(1);(-centeroriginal(3))/(Oaw(3)-centeroriginal(3))*(Oaw(2)-centeroriginal(2))+centeroriginal(2);0];

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

    q1original = [q1new(1)+Oaw(1);q1new(2)+Oaw(2);q1new(3)+Oaw(3)];
    q2original = [q2new(1)+Oaw(1);q2new(2)+Oaw(2);q2new(3)+Oaw(3)];
    q3original = [q3new(1)+Oaw(1);q3new(2)+Oaw(2);q3new(3)+Oaw(3)];
    q4original = [q4new(1)+Oaw(1);q4new(2)+Oaw(2);q4new(3)+Oaw(3)];
    centeroriginal = [centernew(1)+Oaw(1);centernew(2)+Oaw(2);centernew(3)+Oaw(3)];

    q1W(:,step+1) = [(-q1original(3))/(Oaw(3)-q1original(3))*(Oaw(1)-q1original(1))+q1original(1);(-q1original(3))/(Oaw(3)-q1original(3))*(Oaw(2)-q1original(2))+q1original(2);0];
    q2W(:,step+1) = [(-q2original(3))/(Oaw(3)-q2original(3))*(Oaw(1)-q2original(1))+q2original(1);(-q2original(3))/(Oaw(3)-q2original(3))*(Oaw(2)-q2original(2))+q2original(2);0];
    q3W(:,step+1) = [(-q3original(3))/(Oaw(3)-q3original(3))*(Oaw(1)-q3original(1))+q3original(1);(-q3original(3))/(Oaw(3)-q3original(3))*(Oaw(2)-q3original(2))+q3original(2);0];
    q4W(:,step+1) = [(-q4original(3))/(Oaw(3)-q4original(3))*(Oaw(1)-q4original(1))+q4original(1);(-q4original(3))/(Oaw(3)-q4original(3))*(Oaw(2)-q4original(2))+q4original(2);0];
    centerW(:,step+1) = [(-centeroriginal(3))/(Oaw(3)-centeroriginal(3))*(Oaw(1)-centeroriginal(1))+centeroriginal(1);(-centeroriginal(3))/(Oaw(3)-centeroriginal(3))*(Oaw(2)-centeroriginal(2))+centeroriginal(2);0];
    
    step = step+1;

end

%%%%doing some ploting
plotme(kai,kaiTarget,'Pan','PanTarget',n)
ylabel("Pan")
title("Pan")

plotme(fai,faiTarget,'Tilt','TiltTarget',n)
ylabel("Tilt")
title("Tilt")

plotme(kaidot,faidot,'Pandot','Tiltdot',n)
ylabel("velocity")
title("velocity")

plotme(u1,u2,'u1','u2',n)
ylabel("input")
title("input")


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
    title("Tracking in the case with  No Noise")
    view([1 1 1])
    boolean = isinfov(positionTargetInertial(i,:),kai(i),fai(i));
    if boolean == 0        
        fprintf("Target not in Field-of-View!!")
        patch([q1W(1,i),q2W(1,i),q3W(1,i),q4W(1,i)],[q1W(2,i),q2W(2,i),q3W(2,i),q4W(2,i)],'red','FaceAlpha',.3)
        break
    else
        patch([q1W(1,i),q2W(1,i),q3W(1,i),q4W(1,i)],[q1W(2,i),q2W(2,i),q3W(2,i),q4W(2,i)],'green','FaceAlpha',.3)
    end
    
    plot3([positionTargetWorld(i,1)],[positionTargetWorld(i,2)],[positionTargetWorld(i,3)],'r.','markersize',25)
    plot3([Oaw(1)],[Oaw(2)],[Oaw(3)],'b.','markersize',25)
    plot3([q1W(1,i)],[q1W(2,i)],[q1W(3,i)],'k.','markersize',10)
    plot3([q2W(1,i)],[q2W(2,i)],[q2W(3,i)],'k.','markersize',10)
    plot3([q3W(1,i)],[q3W(2,i)],[q3W(3,i)],'k.','markersize',10)
    plot3([q4W(1,i)],[q4W(2,i)],[q4W(3,i)],'k.','markersize',10)
    
    plot3([q1W(1,i),q2W(1,i)],[q1W(2,i),q2W(2,i)],[q1W(3,i),q2W(3,i)],'g','LineWidth',1)
    plot3([q2W(1,i),q3W(1,i)],[q2W(2,i),q3W(2,i)],[q2W(3,i),q3W(3,i)],'g','LineWidth',1)
    plot3([q3W(1,i),q4W(1,i)],[q3W(2,i),q4W(2,i)],[q3W(3,i),q4W(3,i)],'g','LineWidth',1)
    plot3([q4W(1,i),q1W(1,i)],[q4W(2,i),q1W(2,i)],[q4W(3,i),q1W(3,i)],'g','LineWidth',1)
    
    plot3([Oaw(1),q1W(1,i)],[Oaw(2),q1W(2,i)],[Oaw(3),q1W(3,i)],'k-.','LineWidth',1)
    plot3([Oaw(1),q2W(1,i)],[Oaw(2),q2W(2,i)],[Oaw(3),q2W(3,i)],'k-.','LineWidth',1)
    plot3([Oaw(1),q3W(1,i)],[Oaw(2),q3W(2,i)],[Oaw(3),q3W(3,i)],'k-.','LineWidth',1)
    plot3([Oaw(1),q4W(1,i)],[Oaw(2),q4W(2,i)],[Oaw(3),q4W(3,i)],'k-.','LineWidth',1)
    plot3([Oaw(1),centerW(1,i)],[Oaw(2),centerW(2,i)],[Oaw(3),centerW(3,i)],'b','LineWidth',1)
    axis([-5 25 -15 15 0 Oaw(3)])
    %save image to struct M
    M(i) = getframe;
end

%Do the animation and save as MP4
writerObj=VideoWriter('ProjectWithoutNoise.mp4','MPEG-4');
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
    fprintf("x=%d,y=%d",ximageplane,yimageplane);
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
    Uomega=[zeros(1,240) ones(1,40)*pi/2/4 zeros(1,180) ones(1,40)*pi/2/4 zeros(1,100) ones(1,40)*pi/2/4 zeros(1,60)];
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
    iell=[2 10:10:350];
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