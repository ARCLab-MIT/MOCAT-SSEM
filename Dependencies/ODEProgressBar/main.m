function []=main()


options=odeset('OutputFcn',@odeprog,'Events',@odeabort);
[T,S]=ode45(@odedyn,[0:0.001:10],zeros(4,1),options);



function [dS]=odedyn(t,S);
dS=magic(4)*S+rand(4,1);