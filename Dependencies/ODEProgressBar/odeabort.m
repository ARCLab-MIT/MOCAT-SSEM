function [value,isterminal,direction]=odeabort(t,S,varargin)


%Other Events Set Here...ie:
% value(2)=max(abs(S(1:18)))-pi/2;


%Test to see if 'simstop' box is closed
value(1)=double(ishandle(95));
isterminal=[1];
direction=[0];

