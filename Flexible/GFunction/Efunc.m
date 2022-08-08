function out1 = Efunc(A_1,A_2,A_3,E1,E2,E3,I1,I2,I3,Rho1,Rho2,Rho3,dq1,dq2,dq3,dq4,dq5,dq6,dq7,dq8,dq9,dq10,dq11,dq12,g,l1,l2,l3,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12)
%Efunc
%    OUT1 = Efunc(A_1,A_2,A_3,E1,E2,E3,I1,I2,I3,Rho1,Rho2,Rho3,DQ1,DQ2,DQ3,DQ4,DQ5,DQ6,DQ7,DQ8,DQ9,DQ10,DQ11,DQ12,G,L1,L2,L3,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    07-Aug-2022 23:19:45

t2 = cos(q1);
t3 = cos(q2);
t4 = cos(q3);
t5 = sin(q1);
t6 = sin(q2);
t7 = sin(q3);
t8 = dq1.^2;
t9 = dq2.^2;
t10 = dq3.^2;
t11 = dq5.^2;
t12 = l1.^2;
t13 = l2.^2;
t14 = l2.^3;
t15 = l3.^2;
t16 = q4.^2;
t17 = q5.^2;
t18 = q6.^2;
t19 = q7.^2;
t20 = q8.^2;
t21 = q9.^2;
t22 = q10.^2;
t23 = q11.^2;
t24 = q12.^2;
t25 = 1.0./l1;
t26 = 1.0./l2;
t27 = 1.0./l3;
t28 = -q2;
t29 = q1+t28;
t30 = cos(t29);
t31 = sin(t29);
et1 = l2.*t11+(t9.*t14)./3.0+(dq8.^2.*l2)./3.0+(dq7.^2.*t14)./1.05e+2+(dq9.^2.*t14)./1.05e+2+(dq2.*dq7.*t14)./1.5e+1-(dq2.*dq9.*t14)./1.0e+1-(dq7.*dq9.*t14)./7.0e+1+l2.*t8.*t12+l2.*t8.*t17+(l2.*t9.*t20)./3.0+q8.*t9.*t13.*(2.0./3.0)+(t9.*t14.*t19)./1.05e+2+(t9.*t14.*t21)./1.05e+2+dq5.*dq8.*l2.*t30+(dq2.*dq7.*q8.*t13)./1.5e+1-(dq2.*dq8.*q7.*t13)./1.5e+1+(dq2.*dq8.*q9.*t13)./1.0e+1-(dq2.*dq9.*q8.*t13)./1.0e+1+dq2.*dq5.*t13.*t31+(dq5.*dq7.*t13.*t31)./6.0-(dq5.*dq9.*t13.*t31)./6.0+l1.*l2.*q5.*t8.*2.0-(q7.*q9.*t9.*t14)./7.0e+1-dq1.*dq8.*l1.*l2.*t31;
et2 = -dq1.*dq8.*l2.*q5.*t31+dq2.*dq5.*l2.*q8.*t31+dq1.*dq2.*l1.*t13.*t30+(dq1.*dq7.*l1.*t13.*t30)./6.0-(dq1.*dq9.*l1.*t13.*t30)./6.0+dq1.*dq2.*q5.*t13.*t30+(dq1.*dq7.*q5.*t13.*t30)./6.0-(dq2.*dq5.*q7.*t13.*t30)./6.0-(dq1.*dq9.*q5.*t13.*t30)./6.0+(dq2.*dq5.*q9.*t13.*t30)./6.0+dq1.*dq2.*l1.*l2.*q8.*t30+dq1.*dq2.*l2.*q5.*q8.*t30+(dq1.*dq2.*l1.*q7.*t13.*t31)./6.0-(dq1.*dq2.*l1.*q9.*t13.*t31)./6.0+(dq1.*dq2.*q5.*q7.*t13.*t31)./6.0-(dq1.*dq2.*q5.*q9.*t13.*t31)./6.0;
et3 = (A_2.*Rho2.*(et1+et2))./2.0+(A_1.*Rho1.*l1.*(t11.*7.0e+1+t8.*t12.*7.0e+1+t8.*t17.*7.0e+1+dq4.^2.*t12.*2.0+dq6.^2.*t12.*2.0+dq1.*dq4.*t12.*1.4e+1-dq1.*dq6.*t12.*2.1e+1-dq4.*dq6.*t12.*3.0+l1.*q5.*t8.*1.4e+2+t8.*t12.*t16.*2.0+t8.*t12.*t18.*2.0+dq1.*dq4.*l1.*q5.*1.4e+1-dq1.*dq5.*l1.*q4.*1.4e+1+dq1.*dq5.*l1.*q6.*2.1e+1-dq1.*dq6.*l1.*q5.*2.1e+1-q4.*q6.*t8.*t12.*3.0))./4.2e+2+(A_1.*E1.*t17.*t25)./2.0+(A_2.*E2.*t20.*t26)./2.0+(A_3.*E3.*t23.*t27)./2.0;
et4 = (A_3.*Rho3.*l3.*(t10.*t15.*7.0e+1+t10.*t23.*7.0e+1+dq10.^2.*t15.*2.0+dq12.^2.*t15.*2.0+dq11.^2.*7.0e+1+dq3.*dq10.*t15.*1.4e+1-dq3.*dq12.*t15.*2.1e+1-dq10.*dq12.*t15.*3.0+l3.*q11.*t10.*1.4e+2+t10.*t15.*t22.*2.0+t10.*t15.*t24.*2.0+dq3.*dq10.*l3.*q11.*1.4e+1-dq3.*dq11.*l3.*q10.*1.4e+1+dq3.*dq11.*l3.*q12.*2.1e+1-dq3.*dq12.*l3.*q11.*2.1e+1-q10.*q12.*t10.*t15.*3.0))./4.2e+2+E1.*I1.*t25.*(t16+t18+q4.*q6).*2.0+E2.*I2.*t26.*(t19+t21+q7.*q9).*2.0+E3.*I3.*t27.*(t22+t24+q10.*q12).*2.0+(A_1.*Rho1.*g.*l1.*(l1.*t5.*6.0+q5.*t5.*6.0+l1.*q4.*t2-l1.*q6.*t2))./1.2e+1;
et5 = (A_3.*Rho3.*g.*l3.*(l3.*t7.*6.0+q11.*t7.*6.0+l3.*q10.*t4-l3.*q12.*t4))./1.2e+1+(A_2.*Rho2.*g.*l2.*(l1.*t5.*1.2e+1+l2.*t6.*6.0+q5.*t5.*1.2e+1+q8.*t6.*6.0+l2.*q7.*t3-l2.*q9.*t3))./1.2e+1;
out1 = et3+et4+et5;