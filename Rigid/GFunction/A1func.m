function A1 = A1func(l1,l2,l3,q1,q2,q3)
%A1func
%    A1 = A1func(L1,L2,L3,Q1,Q2,Q3)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    07-Aug-2022 19:45:26

A1 = reshape([-l1.*sin(q1),-l2.*sin(q2),l3.*sin(q3),l1.*cos(q1),l2.*cos(q2),-l3.*cos(q3)],[3,2]);
