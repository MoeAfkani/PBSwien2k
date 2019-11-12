%%%
% Split the DOS into two parts gv and gc. Then shift gc to strats from Ec(1)=0;
% and shift and flip gv to start from Ev(1)=0 and increase with energy. Ec
% and Ev are positive and increase from zero. 
% Ev and Ec are in J
% gv and gc are in m^-3.J^-3
% Ref of energy for Ev and Ec are zero; i.e. Ev(1)=0 and Ec(1)=0, Ev is upside down axis
% gv(1)=0 and increases with Ev
% gc(1)=0 and increases with Ec
% Ref of Ef  is Ev(1)
% Ref of Efv is Ev(1), so Efv=Ef
% Ref of Efc is Ec(1), so Efc=-Eg-Ef
% Ed: dopant ionization energy
% Enter [],[] if any of Ev,gv or Ec,gc are not considered
% ND>0 for n type, and ND<0 for p type
% y>0 for electrons and y<0 for holes


function y=FermiDFT(ND,Ev,gv,Ec,gc,Eg,T,Ed,Ef_guess,step,resolution)

kB=1.38047e-23;        %J/K

if isempty(resolution)
    resolution=1e-4;
end

if isempty(Ef_guess)
    Ef=0;
else
    Ef=Ef_guess;
end
if isempty(step) || step<resolution
    step=0.5;
end

Edr=Ed/kB/T;
Efr=Ef/kB/T;

nd=density(Ef,Ev,gv,Ec,gc,Eg,T);

if ND>0; f=2; sg=-1; else f=1/2; sg=1; end %dopant ionization, ref: singh p. 268

d(1)=ND/(1+f*exp(sg*(Edr-Efr)))-nd;
% Solve for Fermi energy:
while abs(step)>resolution
    while d(1)<0
        Efr=Efr+step;
        n=density(Efr*kB*T,Ev,gv,Ec,gc,Eg,T);
        d(1)=ND/(1+f*exp(sg*(Edr-Efr)))-n;
    end
    step=step/2;
    while d(1)>0
        Efr=Efr-step;
        n=density(Efr*kB*T,Ev,gv,Ec,gc,Eg,T);
        d(1)=ND/(1+f*exp(sg*(Edr-Efr)))-n;
    end
    step=step/2;
end

y=-Efr*kB*T;

function y=density(Ef,Ev,gv,Ec,gc,Eg,T)

kB=1.38047e-23;        %J/K

Efv=Ef;
Efc=-Eg-Ef;
Efvr=Efv/kB/T;
Efcr=Efc/kB/T;

if isempty(gv)
    yv=0;
else
    xv=linspace(0,max(0,Efvr)+20,1e3);
    gvi=interp1(Ev/kB/T,gv,xv);
    yiv=gvi./(1+exp(xv-Efvr));
    yv=kB*T*trapz(xv,yiv);
end
if isempty(gc)
    yc=0;
else
    
    xc=linspace(0,max(0,Efcr)+20,1e3);
    gci=interp1(Ec/kB/T,gc,xc);
    yic=gci./(1+exp(xc-Efcr));
    yc=kB*T*trapz(xc,yic);
end
y=-yv+yc;











