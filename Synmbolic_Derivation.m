clc
clear all

syms t
syms X Y Z q1 q2 q3 gama1 gama2 b1 b2% q1 = sai , q2 = theta , q3 = phi , b1 = betaY , b2 = betaZ 
syms dX dY dZ dq1 dq2 dq3 dgama1 dgama2  db1 db2% dq1 = dsai , dq2 = dtheta , dq3 = dphi , db1 = dbetaX , db2 = dbetaZ
syms ms mp mf Ixx Iyy Izz Ixy Ixz Iyz a b c h l ks gamas

% Kinetic Energy of Main Body
VG_XYZ = [dX ; dY ; dZ];
RZ = [cos(q1) sin(q1) 0;-sin(q1) cos(q1) 0; 0 0 1];
RY = [cos(q2) 0 -sin(q2);0 1 0;sin(q2) 0 cos(q2)];
RX = [1 0 0;0 cos(q3) sin(q3);0 -sin(q3) cos(q3)];
R = RX*RY*RZ;
ws_xyz = RX*RY*[0;0;dq1] + RX*[0 ; dq2 ; 0] +  [dq3 ; 0 ; 0]; 
Is = [Ixx -Ixy -Ixz; -Ixy Iyy -Iyz; -Ixz -Iyz Izz];
Hs = Is*ws_xyz;

TSat = 1/2*ms*VG_XYZ.'*VG_XYZ + 1/2*ws_xyz.'*Hs;

% Kinetic Energy of Panel 1
rA1B1m_G = [b/2;-a/2;0];
rGp1_A1B1m = b/2*[-cos(gama1);-sin(gama1);0];
wp1_xyz = ws_xyz + [0;0;dgama1];
Vp1_xyz = R*VG_XYZ + cross(ws_xyz,rA1B1m_G) + cross(wp1_xyz,rGp1_A1B1m);
Ip = 1/12*mp*[c^2 0 0;0 c^2+b^2 0; 0 0 b^2];
R1 = [cos(gama1) sin(gama1) 0;-sin(gama1) cos(gama1) 0; 0 0 1];
wp1_xyzp1 = R1*wp1_xyz;
Hp1 = Ip*wp1_xyzp1;

Tp1 = 1/2*mp*Vp1_xyz.'*Vp1_xyz + 1/2*wp1_xyzp1.'*Hp1;

% Kinetic Energy of Panel 2
rA2B2m_G = [b/2;a/2;0];
rp2_A2B2m = b/2*[-cos(gama2);sin(gama2);0];
wp2_xyz = ws_xyz + [0;0;-dgama2];
Vp2_xyz = R*VG_XYZ + cross(ws_xyz,rA2B2m_G) + cross(wp2_xyz,rp2_A2B2m);
R2 = [cos(gama2) -sin(gama2) 0;sin(gama2) cos(gama2) 0; 0 0 1];
wp2_xyzp2 = R2*wp2_xyz;
Hp2 = Ip*wp2_xyzp2;

Tp2 = 1/2*mp*Vp2_xyz.'*Vp2_xyz + 1/2*wp2_xyzp2.'*Hp2;

% Kinetic Energy of fule
Rz2 = [cos(b2) sin(b2) 0;-sin(b2) cos(b2) 0; 0 0 1];
Rx2 = [1 0 0; 0 cos(b1) sin(b1); 0 -sin(b1) cos(b1)];
rHF_G = [0;-h;0];
rF_HF = (Rz2*Rx2).'*[0;l;0];
wF_xyzF = Rz2*Rx2*ws_xyz + Rz2*[0;0;db2] + [db1;0;0];
wF_xyz = Rx2.'*Rz2.'*wF_xyzF;
VF_xyz =  R*VG_XYZ + cross(ws_xyz,rHF_G) + cross(wF_xyz,rF_HF);

TF = 1/2*mf*VF_xyz.'*VF_xyz;

% Total Kinetic Energy
T = TSat + Tp1 + Tp2 + TF;

% Potential Energy
V = 1/2*ks*(gama1-gamas)^2 + 1/2*ks*(gama2-gamas)^2;

% Lagrangian
L = T - V;

q = [X ; Y ; Z ; q1 ; q2 ; q3 ; gama1 ; gama2 ; b1 ; b2];
dq = [dX ; dY ; dZ ; dq1 ; dq2 ; dq3 ; dgama1 ; dgama2 ; db1 ; db2];   

dL_dq = jacobian(L,q);
dL_ddq = jacobian(L,dq);

M = jacobian(dL_ddq,dq);
N = jacobian(dL_ddq,q);
A = diff(dL_ddq.',t);
B = N*dq + A - dL_dq.';

M = simplify(M);
B = simplify(B);

% Energy 
E = simplify(T + V);

% Linear Momentum
P = simplify(ms*VG_XYZ + R.'*mp*Vp1_xyz + R.'*mp*Vp2_xyz + R.'*mf*VF_xyz);

% Angular Momentum
rG_O = [X;Y;Z];
Hsat_O = R.'*Hs + cross(rG_O,ms*VG_XYZ);

rP1_O = R.'*(rA1B1m_G + rGp1_A1B1m) + rG_O;
HP1_O = R.'*R1.'*Hp1 + cross(rP1_O,mp*R.'*Vp1_xyz);
 
rP2_O = R.'*(rA2B2m_G + rp2_A2B2m) + rG_O;
HP2_O = R.'*R2.'*Hp2 + cross(rP2_O,mp*R.'*Vp2_xyz);

rF_O = R.'*(rHF_G + rF_HF) + rG_O;
HF_O = cross(rF_O,mf*R.'*VF_xyz);

HTot_O = Hsat_O + HP1_O + HP2_O + HF_O;

