clc
clear all
close all
%%
% q1=theta1   q2=theta2 

addpath('GFunction')

syms t x
syms q1 q2 q3
syms dq1 dq2 dq3
syms l1 l2 l3 l4 m1 m2 m3 g 
tor1=0;


%% Link 1
%% Link 1
Q1=[cos(q1) -sin(q1);sin(q1) cos(q1)];
r1(x)=Q1*[x;0];
v1(x)=diff(r1,q1)*dq1;
r1(l1);
v1(l1);
T1=0.5*m1/l1*int(v1.'*v1,x,0,l1);
V1=m1*g/l1*int(r1.'*[0;1],x,0,l1);

%% Link 2
Q2=[cos(q2) -sin(q2);sin(q2) cos(q2)];
r2(x)=r1(l1)+Q2*([x;0]);
v2(x)=diff(r2,q1)*dq1+diff(r2,q2)*dq2;
T2=0.5*m2/l2*int(v2.'*v2,x,0,l2);
V2=m2*g/l2*int(r2.'*[0;1],x,0,l2);

%% Link 3
Q3=[cos(q3) -sin(q3);sin(q3) cos(q3)];
r3(x)=[l4;0]+Q3*([x;0]);
v3(x)=diff(r3,q1)*dq1+diff(r3,q2)*dq2+diff(r3,q3)*dq3;
T3=0.5*m3/l3*int(v3.'*v3,x,0,l3);
V3=m3*g/l3*int(r3.'*[0;1],x,0,l3);


%% Energy
T=T1+T2+T3;
V=V1+V2+V3;
E=T+V;

q=[q1 q2 q3].';
dq=[dq1 dq2 dq3].';

Mass_matrix=jacobian(jacobian(T,dq),dq);
H_matrix=(jacobian(jacobian(T,dq),q))*dq-jacobian(T,q).'+jacobian(V,q).';
M1=[tor1;0];

A2=r2(l2)-r3(l3);
% A=diff(A,q1)*dq1+diff(A,q2)*dq2+diff(A,q3)*dq3+diff(A,q4)*dq4+diff(A,q5)*dq5;
%  
A=jacobian(A2,q)*dq;

A1=jacobian(A,dq).';
% A1=jacobian(A,q).';
AE=simplify(A2.'*A2);


A1dot=A1;
for n=1:size(A1,1)
       A1dot(n,:)= (jacobian(A1(n,:),q)*dq).';
end

M=simplify(Mass_matrix);
H_matrix=simplify(H_matrix);
A1=simplify(A1);
A1dot=simplify(A1dot);
%%

Mfunc=matlabFunction(M,'File','GFunction\Mfunc');
Hfunc=matlabFunction(H_matrix,'File','GFunction\Hfunc');
A1func=matlabFunction(A1,'File','GFunction\A1func');
A1dotfunc=matlabFunction(A1dot,'File','GFunction\A1dotfunc');
Efunc=matlabFunction(simplify(E),'File','GFunction\Efunc');
AEfunc=matlabFunction((AE),'File','GFunction\AEfunc');

