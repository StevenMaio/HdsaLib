clear
clc
close all

theta_bar = [1.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0];
theta_star = [1.1 ; 0.1 ; 0.1 ; 0.1 ; 0.1 ; 0.1 ; 0.1 ; 0.1 ; 0.1];
theta_true = [1.25 ; 0.25 ; 0.25 ; 0.25 ; 0.25 ; 0.25 ; 0.25 ; 0.25 ; 0.25];

writematrix(theta_bar,'theta_bar.txt','Delimiter',' ')
writematrix(theta_star,'theta_star.txt','Delimiter',' ')
writematrix(theta_true,'theta_true.txt','Delimiter',' ')