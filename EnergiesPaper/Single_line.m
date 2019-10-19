%%
%
%decompose a single line from only one side  
%
%
function [Omiga_RPL,Omiga_XPL]=Single_line(T0,IL,U,Z)%,Omiga_RM,Omiga_XM,phi,Omiga_L,RE_D,XE_D,ZE_D

%%Processing the input data
for i=1:size(T0,1)
    Fi(i)=T0(i,2)*pi/180;%The OUTPUT current source phase
    T2(i)=T0(i,1)*exp(1i*Fi(i));%The OUTPUT end current source
end



%%
%Output

%Real and imag part of the line connecting side by side or parallel connection
RE=(real(Z)^2+imag(Z)^2)/real(Z);
XE=(real(Z)^2+imag(Z)^2)/imag(Z);

%The Voltage difference
if (abs(U(1))> abs(U(2)))
    du=U(1)-U(2);
else
    du = U(2) - U(1);
end

%R/X impedance
theta=atan(imag(Z)/real(Z));


% Angle of voltage dif
phi=angle(du);%atan(imag(du)/real(du));

ILR=abs(IL)*cos(theta)*exp(1i*phi);% IL orthogonal decomposition
ILX=abs(IL)*sin(theta)*exp(1i*(phi-pi/2));% IX orthogonal decomposition

%angle dif phi-Fi for projection
angle_dif=phi-Fi;

%Current sets
for i = 1:length(T2)
Omiga_R(i)=abs(T2(i))*cos(phi-Fi(i))*exp(1i*phi);%Real 
Omiga_X(i)=abs(T2(i))*sin(phi-Fi(i))*exp(1i*(phi-pi/2));%Imag 
Omiga_RM(i)=abs(T2(i))*cos(phi-Fi(i));%M means magnitude
Omiga_XM(i)=abs(T2(i))*sin(phi-Fi(i));
%Power of each current in Omiga_R and Omiga_X
%P_Omiga_R(i)=conj(Omiga_R(i)*U(1));%Here P means power on the transmission line
%P_Omiga_X(i)=conj(Omiga_X(i)*U(1));
end

%The positive current set
Omiga_RMP=Omiga_RM;
Omiga_RMP(Omiga_RMP<0)=0;

Omiga_XMP=Omiga_XM;
Omiga_XMP(Omiga_XMP<0)=0;

%Current flowing through the transmission line is the positive part
Omiga_RPL=exp(1i*phi)*abs(IL)*cos(theta)*Omiga_RMP/sum(Omiga_RMP);
Omiga_XPL=exp(1i*(phi-pi/2))*abs(IL)*sin(theta)*Omiga_XMP/sum(Omiga_XMP);

%Decomposition of RE and XE
RE_D=RE*sum(Omiga_RMP)./Omiga_RMP;
XE_D=XE*sum(Omiga_XMP)./Omiga_XMP;

%Parallel connection of RE_D and XE_D
ZE_D=1./(1./RE_D+1./(1i*XE_D));

%Current flowing through ZE_D

Omiga_L = Omiga_RPL + Omiga_XPL;

