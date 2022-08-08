function dZ=AM_dynamics(t,Z)

dZ=Z;


global m1 m2 m3 l1 l2 l3 l4 g
 

q=Z(1:3); dq=Z(4:6);

q1=q(1); q2=q(2); q3=q(3); 

% dq(1)=125;

dq1=dq(1); 
dq2=dq(2); 
dq3=dq(3); 


M = Mfunc(l1,l2,l3,m1,m2,m3,q1,q2);
H = Hfunc(dq1,dq2,g,l1,l2,l3,m1,m2,m3,q1,q2,q3);
A1 = A1func(l1,l2,l3,q1,q2,q3).';
A1dot = A1dotfunc(dq1,dq2,dq3,l1,l2,l3,q1,q2,q3).';


% 
% 
if t<=0.7
%     TO=0;
    TO=1*(1-exp(-t/0.167));
elseif t>0.7
     TO=0.1*(1-exp(-t/0.167));
end

% size(H)
FTT=[TO;0;0];
F=FTT-H;


M_AM=[M -A1.';-A1 zeros(2)];
F_AM=[F;A1dot*dq];

X=pinv(M_AM)*F_AM;

dZ=[dq;X(1:3)];
end