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


function y=SeebeckDFT(Ef,Ev,gv,Ec,gc,Eg,T,r)

kB=1.38047e-23;        %J/K
e=1.602e-19;

Efv=-Ef;
Efc=Ef-Eg;
Efvr=Efv/kB/T;
Efcr=Efc/kB/T;

if isempty(gv)
    yv=0;
else
    xv=linspace(0,max(0,Efvr)+20,1e3);
    gvi=interp1(Ev/kB/T,gv,xv);
    yivt=gvi.*xv.^(1+r).*(xv-Efvr).*exp(xv-Efvr)./(1+exp(xv-Efvr)).^2;
    yivb=gvi.*xv.^(1+r).*exp(xv-Efvr)./(1+exp(xv-Efvr)).^2;
    yvt=trapz(xv,yivt);
    yvb=trapz(xv,yivb);
    yv=kB/e*yvt/yvb;
    yiv=gvi./(1+exp(xv-Efvr));
    p=kB*T*trapz(xv,yiv);
end
if isempty(gc)
    yc=0;
else
    xc=linspace(0,max(0,Efcr)+20,1e3);
    gci=interp1(Ec/kB/T,gc,xc);
    yict=gvi.*xc.^(1+r).*(xc-Efcr).*exp(xc-Efcr)./(1+exp(xc-Efcr)).^2;
    yicb=gci.*xc.^(1+r).*exp(xc-Efcr)./(1+exp(xc-Efcr)).^2;
    yct=trapz(xc,yict);
    ycb=trapz(xc,yicb);
    yc=kB/e*yct/ycb;
    yic=gci./(1+exp(xc-Efcr));
    n=kB*T*trapz(xc,yic);    
end

y=(p*yv-n*yc)/(n+p);











