function AE = AEfunc(l1,l2,l3,l4,q1,q2,q3)
%AEfunc
%    AE = AEfunc(L1,L2,L3,L4,Q1,Q2,Q3)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    07-Aug-2022 19:45:28

AE = (l4-l1.*cos(q1)-l2.*cos(q2)+l3.*cos(q3)).^2+(l1.*sin(q1)+l2.*sin(q2)-l3.*sin(q3)).^2;
