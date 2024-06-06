clear
close all
clc

sabl_evecs = load('Sabl_Output.mat','evecs').evecs;
sabl_evals = load('Sabl_Output.mat','evals').evals;

evals = load('evals.txt')';

evecs = 0*sabl_evecs;
for k = 1:size(evecs,2)
    evecs(:,k) = load(['Evec/Vector_',num2str(k),'.txt']);
end

diff = [norm(evals-sabl_evals);norm(evecs-sabl_evecs)];

figure,
plot(evals,'o')
title('HdsaLib')

figure,
plot(sabl_evals,'o')
title('Sabl')

pause()

for k = 1:size(evecs,2)
    figure,
    plot(linspace(0,1,50),evecs(:,k))
    title('HdsaLib')

    figure,
    plot(linspace(0,1,50),sabl_evecs(:,k))
    title('Sabl')

    pause()
end