f= [1,14,29]
x= [1,2,3] 
y= [0.1,02,0.005] 
z= [1,2,3] 
Yprimenum = diff(f,2)./diff(y,2);
disp( Yprimenum) ;

[x,y] = meshgrid(-4:4,-3:3);
U = x.*x+y.*y*2  
del2(U)

del