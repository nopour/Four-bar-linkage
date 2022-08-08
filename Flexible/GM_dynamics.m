function dZ=GM_dynamics(t,Z)

dZ=Z;


global  m1 m2 m3 l1 l2 l3 l4 g A_1 A_2 A_3 E1 E2 E3 I1 I2 I3 Rho1 Rho2 Rho3
 

q=Z(1:12); dq=Z(13:24);

q1=q(1); q2=q(2); q3=q(3); 
q4=q(4); q5=q(5); q6=q(6); 
q7=q(7); q8=q(8); q9=q(9); 
q10=q(10); q11=q(11); q12=q(12); 


dq1=dq(1); 
dq2=dq(2); 
dq3=dq(3); 
dq4=dq(4); 
dq5=dq(5); 
dq6=dq(6); 
dq7=dq(7); 
dq8=dq(8); 
dq9=dq(9); 
dq10=dq(10); 
dq11=dq(11); 
dq12=dq(12); 


M =  Mfunc(A_1,A_2,A_3,Rho1,Rho2,Rho3,l1,l2,l3,q1,q2,q4,q5,q6,q7,q8,q9,q10,q11,q12);
H =  Hfunc(A_1,A_2,A_3,E1,E2,E3,I1,I2,I3,Rho1,Rho2,Rho3,dq1,dq2,dq3,dq4,dq5,dq6,dq7,dq8,dq9,dq10,dq11,dq12,g,l1,l2,l3,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12);
A1 = A1func(l1,l2,l3,q1,q2,q3,q5,q8,q11).';
A1dot = A1dotfunc(dq1,dq2,dq3,dq5,dq8,dq11,l1,l2,l3,q1,q2,q3,q5,q8,q11).';


% 

if t<=0.7
%     TO=0;
    TO=1*(1-exp(-t/0.167));
elseif t>0.7
     TO=0.1*(1-exp(-t/0.167));
end


% size(H)
% size(M)
% size(A1)
% size(A1dot)

FTT=[TO;0;0;TO;0;0;0;0;0;0;0;0];
F=FTT-H;


F_G=F-A1'*(pinv(A1*(pinv(M)*A1'))*(A1dot*dq+A1*(pinv(M)*F)));

ddq=pinv(M)*F_G;

dZ=[dq;ddq];


end