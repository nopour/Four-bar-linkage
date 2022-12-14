function AE = AEfunc(l1,l2,l3,l4,q1,q2,q3,q5,q8,q11)
%AEfunc
%    AE = AEfunc(L1,L2,L3,L4,Q1,Q2,Q3,Q5,Q8,Q11)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    07-Aug-2022 23:19:46

t2 = l1+q5;
t3 = l2+q8;
t4 = l3+q11;
AE = (l4-t2.*cos(q1)-t3.*cos(q2)+t4.*cos(q3)).^2+(t2.*sin(q1)+t3.*sin(q2)-t4.*sin(q3)).^2;
