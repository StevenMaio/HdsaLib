clear
clc
close all

theta_bar = [1.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0];
theta_star = theta_bar + 0.6;
theta_true = theta_bar + 0.5;

writematrix(theta_bar,'theta_bar.txt','Delimiter',' ')
writematrix(theta_star,'theta_star.txt','Delimiter',' ')
writematrix(theta_true,'theta_true.txt','Delimiter',' ')