% Coordenadas

x1 = 0;
y1 = 0;
z1 = (ma*la + mr*lr)/(ma + mr);

x2 = L2/2 * sin(q2(t))*cos(q1(t));
y2 = L2/2 * sin(q2(t))*sin(q1(t));
z2 = L1 + L2/2*cos(q2(t));

x3 = ( L2 * sin(q2(t)) + L3/2 * sin(q3(t)) ) * cos(q1(t));
y3 = ( L2 * sin(q2(t)) + L3/2 * sin(q3(t)) ) * sin(q1(t));
z3 = L1 + L2/2 * cos(q2(t)) + L3/2 * cos(q3(t));

x4 = x3 + L3/2 * sin(q3(t)) * cos(q1(t)) + ( L4/2 * sin(q4(t)) * cos(q1(t)) );
y4 = y3 + L3/2 * sin(q3(t)) * sin(q1(t)) + ( L4/2 * sin(q4(t)) * sin(q1(t)) );
z4 = z3 + L3/2 * cos(q3(t)) + L4/2 * cos(q4(t));

x5 = x4 + L4/2 * sin(q4(t)) * cos(q1(t)) + ( L5/2 * sin(q5(t)) * cos(q1(t)) );
y5 = y4 + L4/2 * sin(q4(t)) * sin(q1(t)) + ( L5/2 * sin(q5(t)) * sin(q1(t)) );
z5 = z4 + L4/2 * cos(q4(t)) + L5/2 * cos(q5(t));

x6 = x5 + L5/2 * sin(q5(t)) * cos(q1(t)) + ( L6/2 * sin(q6(t)) * cos(q1(t)) );
y6 = y5 + L5/2 * sin(q5(t)) * sin(q1(t)) + ( L6/2 * sin(q6(t)) * sin(q1(t)) );