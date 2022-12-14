function H_matrix = Hfunc1(A_1,A_2,A_3,E1,E2,E3,I1,I2,I3,Rho1,Rho2,Rho3,dq1,dq2,dq3,dq4,dq5,dq6,dq7,dq8,dq9,dq10,dq11,dq12,g,l1,l2,l3,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12)
%Hfunc1
%    H_matrix = Hfunc1(A_1,A_2,A_3,E1,E2,E3,I1,I2,I3,Rho1,Rho2,Rho3,DQ1,DQ2,DQ3,DQ4,DQ5,DQ6,DQ7,DQ8,DQ9,DQ10,DQ11,DQ12,G,L1,L2,L3,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    07-Aug-2022 00:03:20

t2 = cos(q1);
t3 = cos(q2);
t4 = cos(q3);
t5 = sin(q1);
t6 = sin(q2);
t7 = sin(q3);
t8 = dq1.^2;
t9 = dq2.^2;
t10 = dq3.^2;
t11 = l1.^2;
t12 = l1.^3;
t13 = l2.^2;
t14 = l2.^3;
t15 = l3.^2;
t17 = l3.^3;
t18 = 1.0./l1;
t19 = 1.0./l2;
t20 = 1.0./l3;
t21 = -q2;
t16 = t13.^2;
t22 = q1+t21;
t25 = A_2.*Rho2.*g.*t3.*t14.*3.5e+1;
t26 = (A_1.*Rho1.*g.*t2.*t11)./1.2e+1;
t27 = (A_3.*Rho3.*g.*t4.*t15)./1.2e+1;
t23 = cos(t22);
t24 = sin(t22);
t28 = A_2.*Rho2.*dq1.*dq5.*t14.*t23.*7.0e+1;
t29 = A_2.*Rho2.*l1.*t8.*t14.*t24.*3.5e+1;
t30 = A_2.*Rho2.*q5.*t8.*t14.*t24.*3.5e+1;
et1 = A_1.*Rho1.*dq1.*dq5.*t11.*(2.0./3.0)+(A_1.*Rho1.*g.*t2.*t11)./2.0+A_2.*Rho2.*dq1.*dq5.*l1.*l2.*2.0+A_1.*Rho1.*dq1.*dq5.*l1.*q5.*(2.0./3.0)+A_2.*Rho2.*dq1.*dq5.*l2.*q5.*2.0+A_1.*Rho1.*dq1.*dq4.*q4.*t12.*(2.0./1.05e+2)-(A_1.*Rho1.*dq1.*dq4.*q6.*t12)./7.0e+1-(A_1.*Rho1.*dq1.*dq6.*q4.*t12)./7.0e+1+A_1.*Rho1.*dq1.*dq6.*q6.*t12.*(2.0./1.05e+2)+A_2.*Rho2.*g.*l1.*l2.*t2+(A_1.*Rho1.*g.*l1.*q5.*t2)./2.0+A_2.*Rho2.*g.*l2.*q5.*t2-(A_1.*Rho1.*g.*q4.*t5.*t11)./1.2e+1+(A_1.*Rho1.*g.*q6.*t5.*t11)./1.2e+1+(A_2.*Rho2.*l1.*t9.*t13.*t24)./2.0+(A_2.*Rho2.*q5.*t9.*t13.*t24)./2.0+A_2.*Rho2.*dq2.*dq8.*l1.*l2.*t23+A_2.*Rho2.*dq2.*dq8.*l2.*q5.*t23+(A_2.*Rho2.*dq2.*dq7.*l1.*t13.*t24)./6.0-(A_2.*Rho2.*dq2.*dq9.*l1.*t13.*t24)./6.0;
et2 = (A_2.*Rho2.*dq2.*dq7.*q5.*t13.*t24)./6.0-(A_2.*Rho2.*dq2.*dq9.*q5.*t13.*t24)./6.0+(A_2.*Rho2.*l1.*l2.*q8.*t9.*t24)./2.0+(A_2.*Rho2.*l2.*q5.*q8.*t9.*t24)./2.0-(A_2.*Rho2.*l1.*q7.*t9.*t13.*t23)./1.2e+1+(A_2.*Rho2.*l1.*q9.*t9.*t13.*t23)./1.2e+1-(A_2.*Rho2.*q5.*q7.*t9.*t13.*t23)./1.2e+1+(A_2.*Rho2.*q5.*q9.*t9.*t13.*t23)./1.2e+1;
mt1 = [et1+et2;(A_2.*Rho2.*l2.*(dq2.*dq8.*l2.*2.8e+2+dq2.*dq8.*q8.*2.8e+2+g.*l2.*t3.*2.1e+2+g.*q8.*t3.*2.1e+2+dq1.*dq5.*l2.*t23.*4.2e+2+dq2.*dq7.*q7.*t13.*8.0-dq2.*dq7.*q9.*t13.*6.0-dq2.*dq9.*q7.*t13.*6.0+dq2.*dq9.*q9.*t13.*8.0+dq1.*dq5.*q8.*t23.*4.2e+2-g.*l2.*q7.*t6.*3.5e+1+g.*l2.*q9.*t6.*3.5e+1-l1.*l2.*t8.*t24.*2.1e+2-l2.*q5.*t8.*t24.*2.1e+2-l1.*q8.*t8.*t24.*2.1e+2-q5.*q8.*t8.*t24.*2.1e+2+dq1.*dq5.*l2.*q7.*t24.*7.0e+1-dq1.*dq5.*l2.*q9.*t24.*7.0e+1+l1.*l2.*q7.*t8.*t23.*3.5e+1-l1.*l2.*q9.*t8.*t23.*3.5e+1+l2.*q5.*q7.*t8.*t23.*3.5e+1-l2.*q5.*q9.*t8.*t23.*3.5e+1))./4.2e+2];
mt2 = [(A_3.*Rho3.*l3.*(dq3.*dq11.*l3.*2.8e+2+dq3.*dq11.*q11.*2.8e+2+g.*l3.*t4.*2.1e+2+g.*q11.*t4.*2.1e+2+dq3.*dq10.*q10.*t15.*8.0-dq3.*dq10.*q12.*t15.*6.0-dq3.*dq12.*q10.*t15.*6.0+dq3.*dq12.*q12.*t15.*8.0-g.*l3.*q10.*t7.*3.5e+1+g.*l3.*q12.*t7.*3.5e+1))./4.2e+2;t26+(A_1.*Rho1.*((dq1.*dq5.*t11)./1.5e+1-q4.*t8.*t12.*(2.0./1.05e+2)+(q6.*t8.*t12)./7.0e+1))./2.0+E1.*I1.*t18.*(q4.*2.0+q6).*2.0+(A_1.*Rho1.*dq1.*dq5.*t11)./3.0e+1];
mt3 = [t18.*(A_1.*E1.*q5.*-6.0e+1+A_1.*Rho1.*t8.*t12.*2.0e+1+A_1.*Rho1.*dq1.*dq4.*t12.*4.0-A_1.*Rho1.*dq1.*dq6.*t12.*6.0-A_1.*Rho1.*g.*t5.*t11.*3.0e+1+A_2.*Rho2.*l2.*t8.*t11.*6.0e+1+A_1.*Rho1.*q5.*t8.*t11.*2.0e+1-A_2.*Rho2.*g.*l1.*l2.*t5.*6.0e+1+A_2.*Rho2.*l1.*l2.*q5.*t8.*6.0e+1+A_2.*Rho2.*l1.*t9.*t13.*t23.*3.0e+1-A_2.*Rho2.*dq2.*dq8.*l1.*l2.*t24.*6.0e+1+A_2.*Rho2.*dq2.*dq7.*l1.*t13.*t23.*1.0e+1-A_2.*Rho2.*dq2.*dq9.*l1.*t13.*t23.*1.0e+1+A_2.*Rho2.*l1.*l2.*q8.*t9.*t23.*3.0e+1+A_2.*Rho2.*l1.*q7.*t9.*t13.*t24.*5.0-A_2.*Rho2.*l1.*q9.*t9.*t13.*t24.*5.0).*(-1.0./6.0e+1)];
mt4 = [-t26-(A_1.*Rho1.*((dq1.*dq5.*t11)./1.0e+1-(q4.*t8.*t12)./7.0e+1+q6.*t8.*t12.*(2.0./1.05e+2)))./2.0+E1.*I1.*t18.*(q4+q6.*2.0).*2.0-(A_1.*Rho1.*dq1.*dq5.*t11)./2.0e+1;(t19.*(t25+t28-t29-t30+E2.*I2.*q7.*1.68e+3+E2.*I2.*q9.*8.4e+2+A_2.*Rho2.*dq2.*dq8.*t14.*2.8e+1-A_2.*Rho2.*q7.*t9.*t16.*4.0+A_2.*Rho2.*q9.*t9.*t16.*3.0))./4.2e+2;A_2.*t19.*(E2.*q8.*-3.0e+1+Rho2.*t9.*t14.*1.0e+1-Rho2.*g.*t6.*t13.*1.5e+1+Rho2.*q8.*t9.*t13.*1.0e+1+Rho2.*dq2.*dq7.*t14.*2.0-Rho2.*dq2.*dq9.*t14.*3.0+Rho2.*dq1.*dq5.*t13.*t24.*3.0e+1+Rho2.*l1.*t8.*t13.*t23.*1.5e+1+Rho2.*q5.*t8.*t13.*t23.*1.5e+1).*(-1.0./3.0e+1)];
mt5 = [(t19.*(-t25-t28+t29+t30+E2.*I2.*q7.*8.4e+2+E2.*I2.*q9.*1.68e+3-A_2.*Rho2.*dq2.*dq8.*t14.*4.2e+1+A_2.*Rho2.*q7.*t9.*t16.*3.0-A_2.*Rho2.*q9.*t9.*t16.*4.0))./4.2e+2;t27+(A_3.*Rho3.*((dq3.*dq11.*t15)./1.5e+1-q10.*t10.*t17.*(2.0./1.05e+2)+(q12.*t10.*t17)./7.0e+1))./2.0+E3.*I3.*t20.*(q10.*2.0+q12).*2.0+(A_3.*Rho3.*dq3.*dq11.*t15)./3.0e+1;(A_3.*t20.*(E3.*q11.*3.0e+1-Rho3.*t10.*t17.*1.0e+1+Rho3.*g.*t7.*t15.*1.5e+1-Rho3.*q11.*t10.*t15.*1.0e+1-Rho3.*dq3.*dq10.*t17.*2.0+Rho3.*dq3.*dq12.*t17.*3.0))./3.0e+1];
mt6 = [-t27-(A_3.*Rho3.*((dq3.*dq11.*t15)./1.0e+1-(q10.*t10.*t17)./7.0e+1+q12.*t10.*t17.*(2.0./1.05e+2)))./2.0+E3.*I3.*t20.*(q10+q12.*2.0).*2.0-(A_3.*Rho3.*dq3.*dq11.*t15)./2.0e+1];
H_matrix = [mt1;mt2;mt3;mt4;mt5;mt6];
