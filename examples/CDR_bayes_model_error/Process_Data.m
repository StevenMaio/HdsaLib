clear
close all
clc

data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state_lofi = data_obj.data;

data_obj = importdata('control.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;

data_obj = importdata('hifi_state.txt', ' ', 2);  %% we need to skip the first two lines
state_hifi = data_obj.data;

writematrix(control,'control_read.txt')
writematrix(control,'Z.txt')
writematrix(state_hifi-state_lofi,'Y.txt')
writematrix(state_lofi,'state_read.txt')