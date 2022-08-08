function A1dot = A1dotfunc(dq1,dq2,dq3,l1,l2,l3,q1,q2,q3)
%A1dotfunc
%    A1dot = A1dotfunc(DQ1,DQ2,DQ3,L1,L2,L3,Q1,Q2,Q3)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    07-Aug-2022 19:45:26

A1dot = reshape([-dq1.*l1.*cos(q1),-dq2.*l2.*cos(q2),dq3.*l3.*cos(q3),-dq1.*l1.*sin(q1),-dq2.*l2.*sin(q2),dq3.*l3.*sin(q3)],[3,2]);