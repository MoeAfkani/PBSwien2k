%==========================================================================
%                              Initialization
%==========================================================================

versus = 0;  %0 is vs T, 1 is vs ND

if versus == 0;  %vs T
    Xlabel='T (K)';
    plot_func1=str2func('semilogy');
    plot_func2=str2func('plot');
    TT=linspace(305,600,30); 
    T=TT(j)  %K
    Xvar=TT;
%     ND=-5.7e19*1e6; %m^-3
    ND=-2.9e20*1e6; %m^-3 %MnTeLi0.003
elseif versus == 1  %vs doping
    Xlabel='N (cm^-^3)';
    plot_func1=str2func('loglog');
    plot_func2=str2func('semilogx');
    NDD=logspace(log10(1.0e16),log10(1e21),40)*1e6;
    ND=NDD(j)  %m^-3
    T=300;
    Xvar=abs(NDD)/1e6;
end

%==========================================================================
%                              Universal constants
%==========================================================================

global kB h hbar e me NA eps0 c

NA=6.022e23;           %Avogadro's number
kB=1.38047e-23;        %J/K
h=6.626e-34;           %Js
hbar=h/2/pi;           %Js
e=1.602e-19;           %C
me=9.1e-31;            %kg
eps0=8.854187817e-12;  %C^2/j/m
c=3e8;                 %m/s

%==========================================================================
%    Material properties
%==========================================================================

Eg=0.8*e;              %j
Ed=sign(ND)*20*e;      %j dopant activation
r=-1/2;                %Scattering exponent
DFTFile='PM_MnTe_DOS.csv';  %E: 1st column in eV, g: 2nd column in states/eV/atom

Natom=32;              %Number of atoms in the cell
Vcell=(828.7159e-30);  %Cell volume

%==========================================================================
% Read DFT DOS data for MnTeLi 
%==========================================================================

[data]=textread(DFTFile,'','delimiter',',');
E=data(:,1)'*e;  %j
g=data(:,2)'*Natom/Vcell/e; %states/j/m^3 

VBpoint=641;   %Index to VB edge in E vector
CBpoint=642;   %Index to CB edge in E vector

Ev=E(VBpoint)-E(VBpoint:-1:1);
gv=g(VBpoint:-1:1);
Ec=E(CBpoint:end)-E(CBpoint);
gc=g(CBpoint:end);

function y=mv(E,gv)

y= (  gv/(E) * (h^3/4/pi) )^(2/3) / 2
end





