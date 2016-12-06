function F = id_sub(x1)

global M T in Exp Ana sector AA BB CC DD 

C = [x1(5),-x1(5),0,0,0 ; -x1(5),x1(5)+x1(4),-x1(4),0,0 ; 0,-x1(4),x1(4)+x1(3),-x1(3),0 ; 0,0,-x1(3),x1(3)+x1(2),-x1(2) ; 0,0,0,-x1(2),x1(2)+x1(1)];
K = [x1(10),-x1(10),0,0,0 ; -x1(10),x1(10)+x1(9),-x1(9),0,0 ; 0,-x1(9),x1(9)+x1(8),-x1(8),0 ; 0,0,-x1(8),x1(8)+x1(7),-x1(7) ; 0,0,0,-x1(7),x1(7)+x1(6)];

AA= [zeros(5), eye(5) ; -inv(M)*K, -inv(M)*C];
BB= [zeros(5,1) ; inv(M)*[0;1;0;0;0]];
CC= [-inv(M)*K, -inv(M)*C];
DD= inv(M)*[0;1;0;0;0];

Ana = lsim(AA,BB,CC,DD,in,T);

err = Exp(sector,:)-Ana(sector,:);

err2= err.^2;

f   = sum(err2);
F   = sum(f);