%The function will calculet the effective mass of electrones and holes.
%EMg(g,E,MaxE) = (((((g)).^2)./E)*(h^3/4/pi)^2/(8))^(1/3))/me
%   g : array of density of E states 
%   E : array of energy
%   MaxE : maximum energy to be considerd
function m = EMg(g,E,MaxE)
%==========================================================================
%                              Universal constants
%==========================================================================
global h e me
h=6.626e-34;           %Js
me=9.1e-31;            %kg
e=1.602e-19;           %C
%==========================================================================

if isempty(MaxE)
    MaxE=0.5; %ev
end

initP=2;
finalP=2;

while E(finalP)/e < MaxE
    finalP = finalP+1;
end

m=(( (((g(initP:finalP)).^2)./E(initP:finalP))*(h^3/4/pi)^2/(8)  ).^(1/3)) /me  ;

end

