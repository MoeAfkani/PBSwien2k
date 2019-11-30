clear all; clc, tic

% x coordinate is in eV; Y coordinate is in states/eV/atom
% MnTe, 4 atoms, 99.8776 angstrom^3
% Li doped MnTe, 32 atoms,  828.7159 angstrom^3
%----

MATERIALDIR ='MnTe';

%===================================
%   Select the material:
%===================================

% MATERIAL = 'AFM_MnTe';   %AFM bandstructure
MATERIAL = 'PM_MnTe';      %PM  bandstructure

%===================================

path(path,[pwd,'\',MATERIALDIR]);
input_deck = ['input_deck_', MATERIAL];

%========================================================================== 
% Calculated Fermi energy and Seebeck coefficient given gc, gv, ND, and T
% Ref of Ef is Ev(1); Bipolar effect is included
%==========================================================================
j=1;
eval(input_deck);
for j=1:size(Xvar,2)
    eval(input_deck);
    Efermi(j)=FermiDFT(ND,Ev,gv,Ec,gc,Eg,T,Ed,[],[],[]);  %j
    S(j)=SeebeckDFT(Efermi(j),Ev,gv,Ec,gc,Eg,T,r);
end

%==========================================================================
%       Plots
%==========================================================================

FontSize = 20;   
LineWidth = 2;  
MarkerSize = 14;

set(0, 'DefaultTextFontSize', FontSize);
set(0, 'DefaultAxesFontSize', FontSize);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesLineWidth', LineWidth);
set(0, 'DefaultAxesTickLength', [0.015 0.020]);
set(0, 'DefaultAxesBox', 'on');
set(0, 'DefaultLineLineWidth', LineWidth);
set(0, 'DefaultLineMarkerSize', 5);
set(0, 'DefaultPatchLineWidth', LineWidth);
set(0, 'defaultLineMarkerSize', MarkerSize);

%== DOS: g, extracted gc from g, extracted gv from g

figure(1)
% subplot(1,3,1)
% plot(E/e,g/1e6*e); 
% plot(E(VBpoint:-1:1)/e,gv/1e6*e,'r'); hold on
% plot(E(CBpoint:end)/e,gc/1e6*e,'b'); hold off
area(E(VBpoint:-1:1)/e,gv/1e6*e,'FaceColor', [1 0.4 0.6]); hold on
area((E(CBpoint:end)+Eg)/e,gc/1e6*e,'FaceColor',[0 1 1]); %Shifted by Eg to fix the bandgap
xlabel('E(eV)')
ylabel('g (eV^-^1cm^-^3)')
hold off

% subplot(1,3,2)
% plot(Ev/e,gv/1e6*e,'r')
% xlabel('E(eV)')
% ylabel('g_v (eV^-^1cm^-^3)')
% 
% subplot(1,3,3)
% plot(Ec/e,gc/1e6*e,'b')
% xlabel('E(eV)')
% ylabel('g_c (eV^-^1cm^-^3)')

%== S and Ef

figure(2)
plot_func2(Xvar,Efermi/e);
xlabel(Xlabel)
ylabel('E_f-E_v (eV)')
grid on

figure(3)
plot_func2(Xvar,S*1e6);
xlabel(Xlabel)
ylabel('S(\muV/K)')
grid on



% S*1e6
%==========================================================================

toc

%%
z= EMg(gc,Ec,[]);

EkMatrix=readEk('Mn7CrTe/Mn7Cr1Te-scf.energydn_1');

FER=0.6635707483;
Ek=[];
rowNum = 1;
while rowNum<length(EkMatrix)
    if  isnan (EkMatrix(rowNum, 4))
        for inLine=rowNum:rowNum+EkMatrix(rowNum,6)
            tempE=[];
            if abs(EkMatrix(inLine,2)-FER ) < 0.5
                tempE=[tempE,EkMatrix(inLine,2)-FER];
                disp(tempE);
            end
        end
        rowNum=rowNum+EkMatrix(rowNum,6)+1;
    end
    Ek=[Ek;EkMatrix(rowNum, 1:3) , tempE];
end
        
        