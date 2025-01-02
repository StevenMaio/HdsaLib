clear
clc
close all

theta_bar = [1.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0];
theta_star = [1.3 ; 0.3 ; 0.3 ; 0.3 ; 0.3 ; 0.3 ; 0.3 ; 0.3 ; 0.3];

writematrix(theta_bar,'theta_bar.txt','Delimiter',' ')
writematrix(theta_star,'theta_star.txt','Delimiter',' ')