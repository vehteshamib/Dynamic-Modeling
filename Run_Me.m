clc, clear all, close all

global ms mp mf Ixx Iyy Izz Ixy Ixz Iyz a b c h l ks gamas

ms = 2000; mp = 20; mf = 20; Ixx = 1000; Iyy = 1200; Izz = 1400; Ixy = 20; Ixz = -70; Iyz = 40;
a = 2; b = 1; c = 1.4; h=0.2; l=0.6;
ks = 20; gamas = pi/2;
tic
z0 = [zeros(10,1);1000;2000;-1000;0.5;0.3;0.7;0;0;0.5;0.4];
options = odeset('maxstep',0.1);
[t,z] = ode45(@Sat_Dynamics,[0:0.01:10],z0,options);
toc
X = z(:,1); Y = z(:,2); Z = z(:,3); q1 = z(:,4); q2 = z(:,5); q3 = z(:,6);
gama1 = z(:,7); gama2 = z(:,8); b1 = z(:,9); b2 = z(:,10);
dX = z(:,11); dY = z(:,12); dZ = z(:,13); dq1 = z(:,14); dq2 = z(:,15); dq3 = z(:,16);
dgama1 = z(:,17); dgama2 = z(:,18); db1 = z(:,19); db2 = z(:,20);

checking_values = check(t,z);

E_error =checking_values(:,1);
PX_error=checking_values(:,2);
PY_error=checking_values(:,3);
PZ_error=checking_values(:,4);
HX_error=checking_values(:,5);
HY_error=checking_values(:,6);
HZ_error=checking_values(:,7);

figure
hold on
box on
plot(t,X,'r')
plot(t,Y,'g')
plot(t,Z,'b')
ylabel('X & Y & Z','fontsize',18)
xlabel('Time (s)','fontsize',18)
legend('X','Y','Z') 

    pause(5)

figure
hold on
box on
plot(t,q1,'r')
plot(t,q2,'g')
plot(t,q3,'b')
ylabel('q_1 & q_2 & q_3','fontsize',18)
xlabel('Time (s)','fontsize',18)
legend('q_1','q_2','q_3') 

    pause(5)

figure
hold on
box on
plot(t,gama1,'r')
plot(t,gama2,'b')
ylabel('gama_1 & gama_2','fontsize',18)
xlabel('Time (s)','fontsize',18)
legend('gama_1','gama_2') 

    pause(5)

figure
hold on
box on
plot(t,b1,'r')
plot(t,b2,'b')
ylabel('beta_1 & beta_2','fontsize',18)
xlabel('Time (s)','fontsize',18)
legend('beta_1','beta_2') 

    pause(5)

figure
hold on
box on
plot(t,E_error,'b')
ylabel('Energy error','fontsize',18)
xlabel('Time (s)','fontsize',18)
legend('Energy error') 

    pause(5)

figure
hold on
box on
plot(t,PX_error,'r')
plot(t,PY_error,'g')
plot(t,PZ_error,'b')
ylabel('P_x & P_y & P_z Errors','fontsize',18)
xlabel('Time (s)','fontsize',18)
legend('P_x Error','P_y Error','P_z Error') 

    pause(5)

figure
hold on
box on
plot(t,HX_error,'r')
plot(t,HY_error,'g')
plot(t,HZ_error,'b')
ylabel('H_x & H_y & H_z Errors','fontsize',18)
xlabel('Time (s)','fontsize',18)
legend('H_x Error','H_y Error','H_z Error') 

    pause(5)

% from this part on, is the modified version of Dr. Nejat's code.
figure
hold on
Len=length(t);
Xcm=t; Ycm=t; Zcm=t;
for i=1:Len
    R = [                          cos(q1(i))*cos(q2(i)),                           cos(q2(i))*sin(q1(i)),        -sin(q2(i))
        cos(q1(i))*sin(q2(i))*sin(q3(i)) - cos(q3(i))*sin(q1(i)), cos(q1(i))*cos(q3(i)) + sin(q1(i))*sin(q2(i))*sin(q3(i)), cos(q2(i))*sin(q3(i))
        sin(q1(i))*sin(q3(i)) + cos(q1(i))*cos(q3(i))*sin(q2(i)), cos(q3(i))*sin(q1(i))*sin(q2(i)) - cos(q1(i))*sin(q3(i)), cos(q2(i))*cos(q3(i))];
    Rz2 = [cos(b2(i)) sin(b2(i)) 0;-sin(b2(i)) cos(b2(i)) 0; 0 0 1];
    Rx2 = [1 0 0; 0 cos(b1(i)) sin(b1(i)); 0 -sin(b1(i)) cos(b1(i))];
    Rt=R';
    P1 = [X(i);Y(i);Z(i)]+Rt*1/2*[-b;-a;c];
    P2 = [X(i);Y(i);Z(i)]+Rt*1/2*[b;-a;c];
    P3 = [X(i);Y(i);Z(i)]+Rt*1/2*[-b;a;c];
    P4 = [X(i);Y(i);Z(i)]+Rt*1/2*[b;a;c];
    P5 = [X(i);Y(i);Z(i)]+Rt*1/2*[b;a;-c];
    P6 = [X(i);Y(i);Z(i)]+Rt*1/2*[b;-a;-c];
    P7 = [X(i);Y(i);Z(i)]+Rt*1/2*[-b;a;-c];
    P8 = [X(i);Y(i);Z(i)]+Rt*1/2*[-b;-a;-c];
    
    Pg1 = [X(i);Y(i);Z(i)]+Rt*(1/2*[b;-a;0]+b/2*[-cos(gama1(i));-sin(gama1(i));0]);
    Pg2 = [X(i);Y(i);Z(i)]+Rt*(1/2*[b;a;0]+b/2*[-cos(gama2(i));sin(gama2(i));0]);
      
    P9  = [X(i);Y(i);Z(i)]+Rt*(1/2*[b;-a;c]+b*[-cos(gama1(i));-sin(gama1(i));0]);
    P10 = [X(i);Y(i);Z(i)]+Rt*(1/2*[b;-a;-c]+b*[-cos(gama1(i));-sin(gama1(i));0]);
    P11 = [X(i);Y(i);Z(i)]+Rt*(1/2*[b;a;c]+b*[-cos(gama2(i));+sin(gama2(i));0]);
    P12 = [X(i);Y(i);Z(i)]+Rt*(1/2*[b;a;-c]+b*[-cos(gama2(i));+sin(gama2(i));0]);
    
    P13 = [X(i);Y(i);Z(i)]+Rt*[0;-h;0];
    P14 = [X(i);Y(i);Z(i)]+Rt*[0;-h;0]+Rt*(Rz2*Rx2)'*[0;l;0];
    
    V12 = [P1,P2]; V13 = [P1,P3]; V24 = [P2,P4]; V34 = [P3,P4]; V45 = [P4,P5]; V56 = [P5,P6]; V57 = [P5,P7]; V78 = [P7,P8]; V18 = [P1,P8]; V37 = [P3,P7]; V26 = [P2,P6]; V86 = [P8,P6];
    V910 = [P9,P10]; V106 = [P10,P6]; V29 = [P2,P9];
    V411= [P4,P11];   V1112= [P11,P12]; V125= [P12,P5];
    V1314=[P13,P14];
    
    plot3(V12(1,:),V12(2,:),V12(3,:),'k-','linewidth',4)
    hold on
    plot3(V13(1,:),V13(2,:),V13(3,:),'k-','linewidth',4)
    plot3(V24(1,:),V24(2,:),V24(3,:),'k-','linewidth',4)
    plot3(V34(1,:),V34(2,:),V34(3,:),'k-','linewidth',4)
    plot3(V45(1,:),V45(2,:),V45(3,:),'k-','linewidth',4)
    plot3(V56(1,:),V56(2,:),V56(3,:),'k-','linewidth',4)
    plot3(V57(1,:),V57(2,:),V57(3,:),'k-','linewidth',4)
    plot3(V78(1,:),V78(2,:),V78(3,:),'k-','linewidth',4)
    plot3(V18(1,:),V18(2,:),V18(3,:),'k-','linewidth',4)
    plot3(V37(1,:),V37(2,:),V37(3,:),'k-','linewidth',4)
    plot3(V26(1,:),V26(2,:),V26(3,:),'k-','linewidth',4)
    plot3(V86(1,:),V86(2,:),V86(3,:),'k-','linewidth',4)
    
    plot3(V910(1,:),V910(2,:),V910(3,:),'r-','linewidth',3)
    plot3(V106(1,:),V106(2,:),V106(3,:),'r-','linewidth',3)
    plot3(V29(1,:),V29(2,:),V29(3,:),'r-','linewidth',3)
    plot3(V26(1,:),V26(2,:),V26(3,:),'r-','linewidth',3)
    plot3(V411(1,:),V411(2,:),V411(3,:),'r-','linewidth',3)
    plot3(V1112(1,:),V1112(2,:),V1112(3,:),'r-','linewidth',3)
    plot3(V125(1,:),V125(2,:),V125(3,:),'r-','linewidth',3)
    plot3(V45(1,:),V45(2,:),V45(3,:),'r-','linewidth',3)   
    
    plot3(V1314(1,:),V1314(2,:),V1314(3,:),'b-','linewidth',3)
    
    axis equal
    axis([X(i)- 2*a  (X(i)+ 2*a) (Y(i)- 2*a) (Y(i)+2*a)  (Z(i)- 2*a) (Z(i)+2*a)])
    str=['Time = ',num2str(t(i))];
    text(X(i),Y(i),str)
    
    hold off
    pause(0.01)
end



