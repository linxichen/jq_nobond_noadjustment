clear all

A = [4 2 9; 1 2 7; 9 3 0];
B = [7 6 8; 2 5 1; 8 4 1];

[AA,BB,Q,Z] = qz(A,B,'real');
[AAS,BBS,QS,ZS] = ordqz(AA,BB,Q,Z,'udo')
Ggamma = AAS\BBS

inv(Ggamma(2:3,2:3))