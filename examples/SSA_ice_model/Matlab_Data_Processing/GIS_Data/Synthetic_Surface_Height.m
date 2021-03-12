clear
clc
close all

n=81;
[X,Y] = meshgrid(linspace(0,1,n));

syms x y

f1 = @(X,Y) 1.5 + .7*exp(- 28*(X-.45).^2 - 32*(Y-.41).^2 ) + .9*exp(-18.6*(X-.58).^2-20.2*(Y-.56).^2);
z1 = f1(x,y);
f2 = @(X,Y) f1(X,Y) - subs(subs(diff(z1,x),x,0),y,Y)*(X-.5*X.^2);
z2 = f2(x,y);
f3 = @(X,Y) f2(X,Y) - subs(subs(diff(z2,x),x,1),y,Y)*(.5*X.^2);
z3 = f3(x,y);
f4 = @(X,Y) f3(X,Y) - subs(subs(diff(z3,y),y,0),x,X)*(Y-.5*Y.^2);
z4 = f4(x,y);
f5 = @(X,Y) f4(X,Y) - subs(subs(diff(z4,y),y,1),x,X)*(.5*Y.^2);
z5 = f5(x,y);

Z = zeros(n,n);
for i = 1:n
   for j = 1:n
      Z(i,j) = subs(subs(z5,x,X(i,j)),y,Y(i,j)); 
   end
end

surf(X,Y,Z)
shading interp
colorbar
view(2)

Z = Z';
fileID = fopen('Synthetic_Surface_Height.txt', 'w');
fprintf(fileID,'%12.8f\n', Z(:));
fclose(fileID);