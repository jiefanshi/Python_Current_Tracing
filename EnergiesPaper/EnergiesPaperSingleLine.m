clc
clear
close all
v=csvread('C:\Users\LEEPS\Desktop\electric_dissection\EnergiesPaperCsv\V.csv');
z=csvread('C:\Users\LEEPS\Desktop\electric_dissection\EnergiesPaperCsv\Z.csv');
t1=csvread('C:\Users\LEEPS\Desktop\electric_dissection\EnergiesPaperCsv\T1.csv');
t2=csvread('C:\Users\LEEPS\Desktop\electric_dissection\EnergiesPaperCsv\T2.csv');
i=1;
m = 4;
n = 3;
while i <= 1000
    V1 = v(2*i-1,:);
    V2 = v(2*i,:);
    V =[V1(1).*exp(1i*V1(2)*pi/180), V2(1).*exp(1i*V2(2)*pi/180)];
    Z = z(i,1)*(z(i,2)+z(i,3)*1i);
    Tp1 = -[t1(m*i,1)+t1(m*i,2)*1j, t1(m*i-1,1)+t1(m*i-1,2)*1j, t1(m*i-2,1)+t1(m*i-2,2)*1j, t1(m*i-3,1)+t1(m*i-3,2)*1j];
    Tp2 = [t2(n*i,1)+t2(n*i,2)*1j, t2(n*i-1,1)+t2(n*i-1,2)*1j, t2(n*i-2,1)+t2(n*i-2,2)*1j];
    
    IN1= conj(Tp1./V(1));
    IN2= conj(Tp2./V(2));
    
    [theta1, rho1] = cart2pol(real(IN1), imag(IN1));
    [theta2, rho2] = cart2pol(real(IN2), imag(IN2));
    theta1 = theta1 * 180/pi;
    theta2 = theta2 * 180/pi;
    T1 = [rho1; theta1]';
    T2 = [rho2; theta2]';
    IL=(V(2) - V(1)) *20.5^2 /Z;
    [Omiga_RPL1,Omiga_XPL1]=Single_line(T1,IL,V,Z);%,Omiga_RM1,Omiga_XM1,phi1,Omiga_L1,RE_D1,XE_D1,ZE_D1
    [Omiga_RPL2,Omiga_XPL2]=Single_line(T2,IL,V,Z);%,Omiga_RM2,Omiga_XM2,phi2,Omiga_L2,RE_D2,XE_D2,ZE_D2
    %[P_out,P_in,Omiga_in_out,OmigaR_in_out,OmigaX_in_out,RE_D_in_out,XE_D_in_out,ZE_D_in_out]=Single_line_dissection(T1,T2,IL,V,Z);
    %Restore the data
    Omiga_Z = Omiga_RPL2+Omiga_XPL2;
    OZ=Omiga_Z(Omiga_Z ~= 0);
    Data(i,:) = [IL,OZ];
    %[processedDataTheta, processedDataRho] = cart2pol(real(Data), imag(Data));
    %processedData = [processedData; processedDataTheta, processedDataRho];
    i=i+1;
end

%%Data processing
%Switching to theta rho format

[processedDataTheta, processedDataRho] = cart2pol(real(Data), imag(Data));
processedData = [processedDataTheta, processedDataRho];
onlyCurrentData = [processedData(:,1),processedData(:,4)];
lable = ones(1000,1);
lable(1:500) = -1;
processedData = [processedData,lable];
onlyCurrentData = [onlyCurrentData,lable];
% 
%headersProcessedData = {'Line_M','WF_M','Ex_M','Line_A','WF_A','Ex_A','Lable'};
%headersOnlyCurrentData = {'Line_M','Line_A','Lable'};
csvwrite('C:\Users\LEEPS\Desktop\electric_dissection\EnergiesPaperCsv\processedData.csv',processedData);
csvwrite('C:\Users\LEEPS\Desktop\electric_dissection\EnergiesPaperCsv\onlyCurrentData.csv',onlyCurrentData);