clear all
clc
[CRMotored,~,~,~] = convertTDMS(false,'20181016_SI_CRsweep_01.tdms');
[CRData3,~,~,~] = convertTDMS(false,'20181016_SI_CRsweep_04.tdms');%10
 
%% Preliminary
clearvars -except CRMotored CRData3
clc
l=10*0.0254;
B=3.25*0.0254;
L=4.5*0.0254;
a=2.25*0.0254;
Vd=pi/4*B^2*L;
 
CRUMP=CRMotored.Data.MeasuredData(3648).Data;
CRUL3=CRData3.Data.MeasuredData(3648).Data; 
 
d=0; 
for j=1:200; %200 cycles
for i=1:3600; %span the entire length of the vector
    d=d+1;
    CRSM(i,j)=CRUMP(d); 
    CRS3(i,j)=CRUL3(d);
end
end
 
for i=1:3600
   MCRSM(i)=mean(CRSM(i,:));
      MCRS3(i)=mean(CRS3(i,:));
end
 
peak=max(MCRSM); %Finding highest pressure value
loc=0; 
for i=1:3600;
    id=MCRSM(i); %compares every M element to its peak value
    if id==peak; %peak value is found
        loc=i %locates index of highest pressure value
        break
    end
end
CA=[0.2:0.2:720];
TLA = -1.4; 
NCA = CA-(CA(loc)-(TLA)); %Professor Lawler 
loc=loc-(TLA)/0.2;
IVC=loc-135/0.2;
 
for i=1:3600
s(i)=a*cos(NCA(i)*pi/180)+(l^2-a^2*sin(NCA(i)*pi/180)^2)^(1/2);
end
 
VcM=Vd/7;
Vc=Vd/10;
 
 
VM=VcM+(pi*B^2/4)*(l+a-s);
V=Vc+(pi*B^2/4)*(l+a-s);
 
 
IP3=mean(CRData3.Data.MeasuredData(3624).Data)*0.01;
IPM=mean(CRMotored.Data.MeasuredData(3624).Data)*0.01;
for i=1:200
    BPM(i)=mean(CRSM(loc-925:loc-875,i)); %900 indices=180 degrees which is BDC
    BP3(i)=mean(CRS3(loc-925:loc-875,i)); %900 indices=180 degrees which is BDC
end
 
 
for i=1:200
    CRSM(:,i)=CRSM(:,i)+(IPM-BPM(i));
    CRS3(:,i)=CRS3(:,i)+(IP3-BP3(i));
end
for i=1:3600
   MCRSM(i)=mean(CRSM(i,:));
      MCRS3(i)=mean(CRS3(i,:));
end
 
S=loc-120/0.2;
E=loc+120/0.2;
BDC=loc-180/0.2;
 
plot(NCA,CRS3,'-k','Linewidth',0.5)
xlabel('Crank Angle (^{o})')
ylabel('Pressure (bar)')
title('Pressure vs Crank Angle')
hold on
plot(NCA,MCRS3,'-r','LineWidth',2)
 
 
SP3=L*mean(CRData3.Data.MeasuredData(3625).Data)/30;
 
d=0;
 
 
for i=S+1:loc-10  %Pressures should be about the same before combustion
    d=d+1;
I3d(d)=MCRS3(i);
end
j=0;
for k=0:0.01:5 %k will range from 0 to 5 to ensure a correct value is chosen
    d=0;
    j=j+1;
for i=S+1:loc-10
    d=d+1;
    Psm3(d)= MCRS3(IVC)*(V(IVC)/V(i))^k;
end
E3(j)=mean(abs((I3d-Psm3)./I3d)); 
end
 
best3=min(E3); %Lowest error will be used to create simulated pressures
 
 
CAsim=NCA(loc-120/0.2:loc+120/0.2);
 
for i=1:length(E3)
    verif3=E3(i);
    if verif3==best3
        verif3=i
        break
    end
end
k3=(verif3-1)*0.01;
 
%simulated pressures
d=0;
    for i=loc-120/0.2:loc+120/0.2
    d=d+1;
    Psim3(d)= MCRS3(IVC)*(V(IVC)/V(i))^k3;
    end
 
 
figure()   
plot(CAsim,Psim3,'-r',CAsim,MCRS3(loc-120/0.2:loc+120/0.2),'-g')
xlabel('Crank Angle (^{o})')
ylabel('Pressure (bar)')
title('Pressure vs Crank Angle for CR=10')
close all
 
 
d=0;
C1=2.28;
C2=0.00324;
Tref3=mean(CRData3.Data.MeasuredData(3617).Data);
R=0.288;
 
SP=L*mean(CRData3.Data.MeasuredData(3625).Data)/30
for j=1:200
    d=0;
    for i=S:E
        d=d+1;
        w3s(d,j)=C1*SP+C2*Vd*Tref3*(CRS3(i,j)-Psim3(d))/(CRS3(BDC,j)*V(BDC));
    end
end
 
     w3sm=MeanPressures(w3s);
 
 
 
figure() %CR=10
    plot(CAsim,w3s,'-g','LineWidth',0.5)
    hold on
    plot(CAsim,w3sm,'-r','LineWidth',2)
    xlabel('Crank Angle (^{o})')
    ylabel('Velocity (m/s)')
    title('Gas Velocity vs Crank Angle for CR=10')
 
 
AFratio3=mean(CRData3.Data.MeasuredData(3630).Data);
mass=mean(CRData3.Data.MeasuredData(3628).Data)*(1+1/AFratio3)*0.001; 
mfuel=mean(CRData3.Data.MeasuredData(3628).Data)/AFratio3*0.001;
fuelflow=mean(CRData3.Data.MeasuredData(3629).Data)*0.001;
mair=mean(CRData3.Data.MeasuredData(3628).Data)*0.001;
Qlhv=42500000;
Supply=Qlhv*fuelflow; %supply is dQ/dt
 
 
V_int=mass*R*Tref3/MCRS3(IVC)*0.01;
 
Vexh=V(IVC)-V_int;
 
 
mass_res=(MCRS3(IVC)*Vexh*100)/(R*mean(CRData3.Data.MeasuredData(3614).Data));
 
 
mass_tot=mass_res+mass;
 
 
d=0;
for j=1:200
    d=0;
for i=S:E
    d=d+1;
Tbulk(d,j)=CRS3(i,j)*V(i)*100/(mass_tot*R);
end 
end
 
for j=1:200
    d=0;
    for i=S:E
        d=d+1;
        h(d,j)=3.26*B^(-0.2)*(CRS3(i,j)*100)^(0.8)*Tbulk(d,j)^(-0.55)*w3s(d,j)^0.8;
    end
end
 
d=0;
for i=1:length(h)
    d=d+1;
    hm(d)=mean(h(i,:));
end
    
figure()    
    plot(CAsim,h,'-g','LineWidth',0.5)
 xlabel('Crank Angle (^{o})')
    ylabel('h (W/(m^{2}-K)')
    title('h vs Crank Angle for CR=10')
    hold on
    plot(CAsim,hm,'-r','LineWidth',2)
           
 
    
A=V/(pi/4*B^2)*pi*B; %h=V/A, A=pi*B*h
Twall=400
 
h=h*1/(1000/60*360); %conversion factors
hm=hm*1/(1000/60*360);
 
 
d=0;
for j=1:200
    d=0;
    for i=S:E
        d=d+1;
       hLoss(d,j)=h(d,j)*A(i)*(Tbulk(d,j)-400);
    end
end
d=0;
hLossm=MeanPressures(hLoss);
 
 
figure()
plot(CAsim,hLoss,'-g','LineWidth',0.5)
xlabel('Crank Angle (^{o})')
ylabel('Heat Loss Rate to Wall (J/CA)')
title('Heat Loss to Wall Rate vs Crank Angle for CR=10')
hold on
plot(CAsim,hLossm,'-r','LineWidth',2)
    
    for j=1:200 %indefinite integral
        d=0;
        hRelease(1,j)=0;
        
        for i=2:length(hLoss)-1
    hRelease(i,j)= hRelease(i-1,j)+0.2*(hLoss(i+1,j)+hLoss(i,j))/2;
        end
    end
hReleasem=MeanPressures(hRelease);
 
    
     figure()
    plot(CAsim(1:1200),hRelease,'-g','LineWidth',0.5)
    xlabel('Crank Angle (^{o})')
    ylabel('Heat Lost (J)')
    title('Heat Lost to Wall for CR=10')
    hold on
    plot(CAsim(1:1200),hReleasem,'-r','LineWidth',2)
       
d=0;
g=1.33; %heat release rates
for j=1:200
    d=0;
  for  i=S:E
    d=d+1;
    Loss(d,j)=(g*100000*CRS3(i,j)*((V(i+1)-V(i-1))/0.4))/(g-1)+(V(i)*100000*((CRS3(i+1,j)-CRS3(i-1,j))/0.4))/(g-1);
  end
end
Lossm=MeanPressures(Loss);
 
figure()
plot(CAsim,Loss,'-g','LineWidth',0.5)
xlabel('Crank Angle (^{o})')
ylabel('Heat Release Rate (J/CA)')
title('Heat Release Rate vs Crank Angle')
hold on
plot(CAsim,Lossm,'-r','LineWidth',2)
 
figure() %Verifying that COV is greatest after combustion
for i=1:length(Loss)
    COVLoss(i)=abs(std(Loss(i,:))/mean(Loss(i,:)));
end
plot(CAsim,COVLoss,'-r')
xlabel('Crank Angle (^{o})')
ylabel('Coefficent of Variation')
title('Coefficient of Variation vs Crank Angle')
 
 
 
for j=1:200 %Heat released 
    Release(1,j)=0;
  for i=2:length(Loss)-1
      Release(i,j)= Release(i-1,j)+0.2*(Loss(i,j)+Loss(i+1,j))/2;
  end 
end
 
Releasem=MeanPressures(Release);
 
TotalHeat=Release+hRelease; %total heat released
TotalHeatm=MeanPressures(TotalHeat);
figure()
plot(CAsim(1:1200),TotalHeat,'-g','LineWidth',0.5)
xlabel('Crank Angle (^{o})')
ylabel('Heat Released (J)')
title('Net Heat Released vs Crank Angle')
hold on
plot(CAsim(1:1200),TotalHeatm,'-r','LineWidth',2)
for j=1:200 %normalize by max to creat MFB Curve
MFB(:,j)=TotalHeat(:,j)/max(TotalHeat(:,j));
end
 
 
MFBm=MeanPressures(MFB);
figure()
plot(CAsim(1:1200),MFB,'-g','LineWidth',0.5)
xlabel('Crank Angle (^{O})')
ylabel('Mass Fraction Burned')
title('Mass Fraction Burned vs Crank Angle')
hold on
plot(CAsim(1:1200),MFBm,'-r','LineWidth',2)
 
w=mean(CRData3.Data.MeasuredData(3625).Data)/60; %omega
T=mean(CRData3.Data.MeasuredData(3626).Data); %Torque
BrakePower=T*2*pi*w;
 
NetWorkMean=0;
for i=1:3599
    NetWorkMean=NetWorkMean+(V(i+1)-V(i))*100000*(MCRS3(i)+MCRS3(i+1))/2;
end
GrossWorkMean=0;
for i=loc-900:loc+900
    GrossWorkMean=GrossWorkMean+(V(i+1)-V(i))*100000*(MCRS3(i)+MCRS3(i+1))/2;
end
 
INPMean=NetWorkMean*w/2;
IGPMean=GrossWorkMean*w/2;
FPMean=IGPMean-BrakePower;
 
INEMean=INPMean/Supply*100;
IGEMean=IGPMean/Supply*100;
BEMean=BrakePower/Supply*100;
MEMean=BrakePower/IGPMean*100;
NFCEMean=NetWorkMean/(Qlhv*mfuel)*100;
 
MeanEfficiencies=[INEMean IGEMean BEMean MEMean NFCEMean]
 
BSFC=mean(CRData3.Data.MeasuredData(3629).Data)*3600*1000/BrakePower;
NSFC=mean(CRData3.Data.MeasuredData(3629).Data)*3600*1000/INPMean;
GSFC=mean(CRData3.Data.MeasuredData(3629).Data)*3600*1000/IGPMean;
 
MeanSFC=[NSFC GSFC BSFC]
 
 
 
 
 
 
 
 
 
 
 
 
 
 
%% Identifying Knock
close all
 
 
for j=1:200 
    for i=1:length(MFB)
            H=MFB(i,j);
        if H>0.5
            z(j)=1/(MFB(i,j)-MFB(i-1,j))*(0.5-MFB(i-1,j))+i-1; %Linear Interpolation
            break
        end
    end
end
 
for i=1:length(z)
CA50(i)=(loc-600+z(i))*0.2-360; %Converting to CA
end
min50=min(CA50)
max50=max(CA50)
for i=1:200
    idB=CA50(i);
    if idB==max50;
        Rb=i; %R is the most retarded case
        break
    end
end
for i=1:200
    idS=CA50(i);
    if idS==min50;
        Ra=i; %Ra is the advanced case
        break
    end
end
figure()
plot(NCA(loc-30/0.2:E),CRS3(loc-30/0.2:E,Ra),'-r',NCA(loc-30/0.2:E),CRS3(loc-30/0.2:E,Rb),'-k')
xlabel('Crank Angle (^{o})')
ylabel('Pressure (bar)')
title('Pressure vs Crank Angle for Most Retarded and Advanced Cases')
legend('Advanced Case','Retarded Case')
 
for j=1:200
    d=0;
for i=S:E
    d=d+1;
    TI(d,j)=Tbulk(1,j)*(CRS3(i,j)/CRS3(S,j))^((g-1)/g); %Isentropic Unburened Temps
end
end
 
for j=1:200
    d=0;
    for i=S:E
        d=d+1;
        tau(d,j)=(1.3*10^(-4))*(CRS3(i,j)*0.986923)^(-1.05)*20.5^(-1.41)*exp(33700/(1.987*TI(d,j)));
    end
end
 
CF=w*360/1000
tau=tau*CF; %convert tau
for j=1:200
    for i=1:length(tau)
        invTau(i,j)=1/(tau(i,j)); %inv tau will be numerically integrated
    end
end
 
for j=1:200
    KO(1,j)=0;
   for i=2:length(tau)-1
    KO(i,j)=KO(i-1,j)+0.2*(invTau(i,j)+invTau(i+1,j))/2; %integrate over CA space
   end
end
 
plot(CAsim(1:1200),KO,'-r')
xlabel('Crank Angle ^({o})')
ylabel('Integral of 1/\tau)')
 
 
 
r=1;
q=1;
for j=1:200
    id=max(KO(:,j)); %Find all knocking and not knocking cases
    if id >= 1
        Knock(r)=j;
        r=r+1;
    else 
        Safe(q)=j
        q=q+1;
    end
end
for j=1:length(Knock)
    for i=1:3600
CRKnock(i,j)=CRS3(i,Knock(j)); %Isolating Knock Cycles
    end
end
for j=1:length(Safe)
    for i=1:3600
    CRSafe(i,j)=CRS3(i,Safe(j)); %Isolating No-Knock Cycles
    end
end
plot(NCA,CRSafe,'-r')
xlabel('Crank Angle (^{o})')
ylabel('Pressure (bar)')
title('Non-Knocking Pressures vs Crank Angle')
xlim([0 35])
ylim([16 35])
 
 
plot(NCA,CRKnock,'-r')
xlabel('Crank Angle (^{o})')
ylabel('Pressure (bar)')
title('Knocking Pressures vs Crank Angle')
xlim([0 35])
ylim([16 35])
 
MCRKnock=MeanPressures(CRKnock);
MCRSafe=MeanPressures(CRSafe);
 
%% Knock and No-Knock Cases
 
% all processes repeated for isolated knock and no knock cases
d=0;
for j=1:163
    d=0;
for i=S:E
    d=d+1;
TbulkK(d,j)=CRKnock(i,j)*V(i)*100/(mass_tot*R);
end 
end
 
d=0;
for j=1:37
    d=0;
for i=S:E
    d=d+1;
TbulkS(d,j)=CRSafe(i,j)*V(i)*100/(mass_tot*R);
end 
end
 
d=0;
C1=2.28;
C2=0.00324;
Tref3=mean(CRData3.Data.MeasuredData(3617).Data);
R=0.288;
 
SP=L*mean(CRData3.Data.MeasuredData(3625).Data)/30
for j=1:163
    d=0;
    for i=S:E
        d=d+1;
        wK(d,j)=C1*SP+C2*Vd*Tref3*(CRKnock(i,j)-Psim3(d))/(CRKnock(BDC,j)*V(BDC));
    end
end
wKm=MeanPressures(wK);
 
for j=1:37
    d=0;
    for i=S:E
        d=d+1;
        wS(d,j)=C1*SP+C2*Vd*Tref3*(CRSafe(i,j)-Psim3(d))/(CRSafe(BDC,j)*V(BDC));
    end
end
wSm=MeanPressures(wS);
 
for j=1:163
    d=0;
    for i=S:E
        d=d+1;
        hK(d,j)=3.26*B^(-0.2)*(CRKnock(i,j)*100)^(0.8)*TbulkK(d,j)^(-0.55)*wK(d,j)^0.8;
    end
end
hKm=MeanPressures(hK);
 
for j=1:37
    d=0;
    for i=S:E
        d=d+1;
        hS(d,j)=3.26*B^(-0.2)*(CRSafe(i,j)*100)^(0.8)*TbulkS(d,j)^(-0.55)*wS(d,j)^0.8;
    end
end
hSm=MeanPressures(hS);
           
 
    
A=V/(pi/4*B^2)*pi*B;
Twall=400
 
hS=hS*1/(1000/60*360);
hSm=hSm*1/(1000/60*360);
 
hK=hK*1/(1000/60*360);
hKm=hKm*1/(1000/60*360);
 
d=0;
for j=1:163
    d=0;
    for i=S:E
        d=d+1;
       hLossK(d,j)=hK(d,j)*A(i)*(TbulkK(d,j)-400);
    end
end
d=0;
hLossKm=MeanPressures(hLossK);
 
d=0;
for j=1:37
    d=0;
    for i=S:E
        d=d+1;
       hLossS(d,j)=hS(d,j)*A(i)*(TbulkS(d,j)-400);
    end
end
d=0;
hLossSm=MeanPressures(hLossS);
 
 
 
 
figure()
plot(CAsim,hLossK,'-g','LineWidth',0.5)
xlabel('Crank Angle (^{o})')
ylabel('Heat Loss Rate to Wall (J/CA)')
title('Heat Loss to Wall Rate vs Crank Angle for CR=10')
hold on
plot(CAsim,hLossKm,'-r','LineWidth',2)
    
    for j=1:163
        d=0;
        hReleaseK(1,j)=0;
        
        for i=2:length(hLossK)-1
    hReleaseK(i,j)= hReleaseK(i-1,j)+0.2*(hLossK(i+1,j)+hLossK(i,j))/2;
        end
    end
hReleaseKm=MeanPressures(hReleaseK);
 
    for j=1:37
        d=0;
        hReleaseS(1,j)=0;
        
        for i=2:length(hLossK)-1
    hReleaseS(i,j)= hReleaseS(i-1,j)+0.2*(hLossS(i+1,j)+hLossS(i,j))/2;
        end
    end
hReleaseSm=MeanPressures(hReleaseS);
 
    
     figure()
    plot(CAsim(1:1200),hReleaseS,'-g','LineWidth',0.5)
    xlabel('Crank Angle (^{o})')
    ylabel('Heat Lost (J)')
    title('Heat Lost to Wall for CR=10')
    hold on
    plot(CAsim(1:1200),hReleaseSm,'-r','LineWidth',2)
    
       
d=0;
g=1.33;
for j=1:163
    d=0;
  for  i=S:E
    d=d+1;
    LossK(d,j)=(g*100000*CRKnock(i,j)*((V(i+1)-V(i-1))/0.4))/(g-1)+(V(i)*100000*((CRKnock(i+1,j)-CRKnock(i-1,j))/0.4))/(g-1);
  end
end
LossKm=MeanPressures(LossK);
 
for j=1:37
    d=0;
  for  i=S:E
    d=d+1;
    LossS(d,j)=(g*100000*CRSafe(i,j)*((V(i+1)-V(i-1))/0.4))/(g-1)+(V(i)*100000*((CRSafe(i+1,j)-CRSafe(i-1,j))/0.4))/(g-1);
  end
end
LossSm=MeanPressures(LossS);
 
figure()
plot(CAsim,LossK,'-g','LineWidth',0.5)
xlabel('Crank Angle (^{o})')
ylabel('Heat Release Rate (J/CA)')
title('Heat Release Rate vs Crank Angle')
hold on
plot(CAsim,LossKm,'-r','LineWidth',2)
 
 
for i=1:length(LossK)
    COVLossK(i)=abs(std(LossK(i,:))/mean(LossK(i,:)));
end
 
for i=1:length(LossS)
    COVLossS(i)=abs(std(LossS(i,:))/mean(LossS(i,:)));
end
figure()
plot(CAsim,COVLossS,'-r')
xlabel('Crank Angle (^{o})')
ylabel('Coefficent of Variation')
title('Coefficient of Variation vs Crank Angle')
 
 
 
for j=1:163
    ReleaseK(1,j)=0;
  for i=2:length(Loss)-1
      ReleaseK(i,j)= ReleaseK(i-1,j)+0.2*(LossK(i,j)+LossK(i+1,j))/2;
  end 
end
 
ReleaseKm=MeanPressures(ReleaseK);
 
for j=1:37
    ReleaseS(1,j)=0;
  for i=2:length(Loss)-1
      ReleaseS(i,j)= ReleaseS(i-1,j)+0.2*(LossS(i,j)+LossS(i+1,j))/2;
  end 
end
 
ReleaseSm=MeanPressures(ReleaseS);
 
TotalHeatK=ReleaseK+hReleaseK;
TotalHeatKm=MeanPressures(TotalHeatK);
TotalHeatS=ReleaseS+hReleaseS;
TotalHeatSm=MeanPressures(TotalHeatS);
 
figure()
plot(CAsim(1:1200),TotalHeatS,'-g','LineWidth',0.5)
xlabel('Crank Angle (^{o})')
ylabel('Heat Released (J)')
title('Net Heat Released vs Crank Angle')
hold on
plot(CAsim(1:1200),TotalHeatSm,'-r','LineWidth',2)
for j=1:163
MFBK(:,j)=TotalHeatK(:,j)/max(TotalHeatK(:,j));
end
 
for j=1:37
MFBS(:,j)=TotalHeatS(:,j)/max(TotalHeatS(:,j));
end
 
 
MFBKm=MeanPressures(MFBK);
MFBSm=MeanPressures(MFBS);
 
figure()
plot(CAsim(1:1200),MFBS,'-g','LineWidth',0.5)
xlabel('Crank Angle (^{O})')
ylabel('Mass Fraction Burned')
title('Mass Fraction Burned vs Crank Angle')
hold on
plot(CAsim(1:1200),MFBSm,'-r','LineWidth',2)
 
T=mean(CRData3.Data.MeasuredData(3626).Data);
BrakePower=T*2*pi*w;
 
NetWorkMeanK=0;
for i=1:3599
    NetWorkMeanK=NetWorkMeanK+(V(i+1)-V(i))*100000*(MCRKnock(i)+MCRKnock(i+1))/2;
end
GrossWorkMeanK=0;
for i=loc-900:loc+900
    GrossWorkMeanK=GrossWorkMeanK+(V(i+1)-V(i))*100000*(MCRKnock(i)+MCRKnock(i+1))/2;
end
 
NetWorkMeanS=0;
for i=1:3599
    NetWorkMeanS=NetWorkMeanS+(V(i+1)-V(i))*100000*(MCRSafe(i)+MCRSafe(i+1))/2;
end
GrossWorkMeanS=0;
 
for i=loc-900:loc+900
    GrossWorkMeanS=GrossWorkMeanS+(V(i+1)-V(i))*100000*(MCRSafe(i)+MCRSafe(i+1))/2;
end
 
INPK=NetWorkMeanK*w/2;
IGPK=GrossWorkMeanK*w/2;
FPK=IGPK-BrakePower;
 
INEK=INPK/Supply*100;
IGEK=IGPK/Supply*100;
BEMean=BrakePower/Supply*100;
MEK=BrakePower/IGPK*100;
NFCEK=NetWorkMeanK/(Qlhv*mfuel)*100;
 
KnockEfficiencies=[INEK IGEK BEMean MEK NFCEK]
 
BSFCK=mean(CRData3.Data.MeasuredData(3629).Data)*3600*1000/BrakePower;
NSFCK=mean(CRData3.Data.MeasuredData(3629).Data)*3600*1000/INPK;
GSFCK=mean(CRData3.Data.MeasuredData(3629).Data)*3600*1000/IGPK;
 
KnockSFC=[NSFCK GSFCK BSFCK]
 
 
INPS=NetWorkMeanS*w/2;
IGPS=GrossWorkMeanS*w/2;
FPS=IGPS-BrakePower;
 
INES=INPS/Supply*100;
IGES=IGPS/Supply*100;
BEMean=BrakePower/Supply*100;
MES=BrakePower/IGPS*100;
NFCES=NetWorkMeanS/(Qlhv*mfuel)*100;
 
NonKnockEfficiencies=[INES IGES BEMean MES NFCES]
 
BSFCS=mean(CRData3.Data.MeasuredData(3629).Data)*3600*1000/BrakePower;
NSFCS=mean(CRData3.Data.MeasuredData(3629).Data)*3600*1000/INPS;
GSFCS=mean(CRData3.Data.MeasuredData(3629).Data)*3600*1000/IGPS;
 
NonKnockSFC=[NSFCS GSFCS BSFCS]
 
Efficiencies=[MeanEfficiencies' KnockEfficiencies' NonKnockEfficiencies'] 
 
SFC=[MeanSFC' KnockSFC' NonKnockSFC']
 
close all
 
for j=1:163
    d=0;
for i=S:E
    d=d+1;
    TIK(d,j)=TbulkK(1,j)*(CRKnock(i,j)/CRKnock(S,j))^((g-1)/g);
end
end
TIKm=MeanPressures(TIK);
 
for j=1:37
    d=0;
for i=S:E
    d=d+1;
    TIS(d,j)=TbulkS(1,j)*(CRSafe(i,j)/CRSafe(S,j))^((g-1)/g);
end
end
 
TIm=MeanPressures(TI);
 
TISm=MeanPressures(TIS);
 
plot(CAsim,TIKm,CAsim,TISm,'-r',CAsim,TIm,'-g')
close all
 
%% Graphs
figure()
p1=plot(NCA,CRS3,'-b','LineWidth',0.3) %PvsCA for all
hold on
p2=plot(NCA,MCRS3,'--r','LineWidth',1.5)
hold on
p3=plot(NCA,MCRSafe,'-k','LineWidth',2)
hold on
p4=plot(NCA,MCRKnock,'-y','LineWidth',2)
hold on
legend([p1(1) p2(1) p3(1) p4(1)],'All cycles','Mean Cycle','No-Knock Cycle','Knock Cycle')
xlabel('Crank Angle (^{o})')
ylabel('Pressure (bar)')
xlim([5 30])
 
figure()
p1=plot(NCA,CRSafe,'-g','LineWidth',0.3) %PvsCA noKnock
hold on
p2=plot(NCA,MCRSafe,'-r','LineWidth',2)
xlabel('Crank Angle (^{o})')
ylabel('Pressure (bar)')
xlim([10 30])
ylim([20 32])
legend([p1(1) p2(1)],'All Cycles','Average Cycle')
 
figure()
p1=plot(NCA,CRKnock,'-g','LineWidth',0.3) %PvsCA Knock
hold on
p2=plot(NCA,MCRKnock,'-r','LineWidth',2)
xlabel('Crank Angle (^{o})')
ylabel('Pressure (bar)')
xlim([10 30])
ylim([20 35])
legend([p1(1) p2(1)],'All Cycles','Average Cycle')
 
 
figure()
loglog(V,MCRS3,'--r') %PV 
hold on
loglog(V,MCRSafe,'-k')
hold on
loglog(V,MCRKnock,'-b')
hold on
legend('Mean Cycle','No-Knock Cycle','Knock Cycle')
xlabel('Volume (m^{3})')
ylabel('Pressure (bar)')
 
figure()
plot(CAsim,Psim3,'-r',CAsim,MCRS3(loc-120/0.2:loc+120/0.2),'-g') %Psim
xlabel('Crank Angle (^{o})')
ylabel('Pressure (bar)')
title('Pressure vs Crank Angle for CR=10')
 
figure()
p1=plot(CAsim,w3s,'-c','LineWidth',0.3) %Gas Velocities
hold on
p2=plot(CAsim,w3sm,'--r','LineWidth',1.5)
hold on
p3=plot(CAsim,wKm,'-k','LineWidth',2)
hold on
p4=plot(CAsim,wSm,'-y','LineWidth',2)
hold on
xlabel('Crank Angle (^{o})')
ylabel('Gas Velocity (m/s)')
legend([p1(1) p2(1) p3(1) p4(1)],'All Cycles','Mean Cycle','Knock Cycle','No-Knock Case')
 
figure()    
Tbulkm=MeanPressures(Tbulk);
TbulkKm=MeanPressures(TbulkK);
TbulkSm=MeanPressures(TbulkS);
 
figure()
p1=plot(CAsim,Tbulk,'-m','LineWidth',0.3) %T bulk
hold on
p2=plot(CAsim,Tbulkm,'--c','LineWidth',1.5)
hold on
p3=plot(CAsim,TbulkKm,'-k','LineWidth',2)
hold on
p4=plot(CAsim,TbulkSm,'-y','LineWidth',2)
hold on
xlabel('Crank Angle (^{o})')
ylabel('Bulk Temperature (K)')
legend([p1(1) p2(1) p3(1) p4(1)],'All Cycles','Mean Cycle','Knock Cycle','No-Knock Case')
 
figure()
p1=plot(CAsim,TbulkK,'-m','LineWidth',0.3) %T bulk K
hold on
p2=plot(CAsim,TbulkKm,'c','LineWidth',1.5)
hold on
xlabel('Crank Angle (^{o})')
ylabel('Bulk Temperature (K)')
legend([p1(1) p2(1)],'All Cycles','Average Cycle')
 
figure()
p1=plot(CAsim,TbulkS,'-m','LineWidth',0.3) %T bulk S
hold on
p2=plot(CAsim,TbulkSm,'c','LineWidth',1.5)
hold on
xlabel('Crank Angle (^{o})')
ylabel('Bulk Temperature (K)')
legend([p1(1) p2(1)],'All Cycles','Average Cycle')
 
 
figure()        
p1=plot(CAsim,h,'-m','LineWidth',0.3) %h
hold on
p2=plot(CAsim,hm,'--c','LineWidth',1.5)
hold on
p3=plot(CAsim,hKm,'-k','LineWidth',2)
hold on
p4=plot(CAsim,hSm,'-y','LineWidth',2)
hold on
xlabel('Crank Angle (^{o})')
ylabel('h (J/(CA-m^{2}-K')
legend([p1(1) p2(1) p3(1) p4(1)],'All Cycles','Mean Cycle','Knock Cycle','No-Knock Case')
 
figure()
p1=plot(CAsim,hK,'-k','LineWidth',0.3) %T bulk S
hold on
p2=plot(CAsim,hKm,'y','LineWidth',1.5)
hold on
xlabel('Crank Angle (^{o})')
ylabel('h (J/(CA-m^{2}-K')
legend([p1(1) p2(1)],'All Cycles','Average Cycle')
xlim([0 30])
 
figure()
p1=plot(CAsim,hS,'-k','LineWidth',0.3) %T bulk S
hold on
p2=plot(CAsim,hSm,'y','LineWidth',1.5)
hold on
xlabel('Crank Angle (^{o})')
ylabel('h (J/(CA-m^{2}-K')
legend([p1(1) p2(1)],'All Cycles','Average Cycle')
xlim([0 30])
 
figure()
p1=plot(CAsim,hLossm,'-r')
hold on
p2=plot(CAsim,hLossKm,'-k')
hold on
p3=plot(CAsim,hLossSm,'-m')
xlabel('Crank Angle (^{o})')
ylabel('Heat Transfer Rate (J/CA)')
legend([p1(1) p2(1) p3(1)],'Mean Cycle','Knock Cycle','No-Knock Cycle')
 
figure()
p1=plot(CAsim(1:1200),hReleasem,'-r')
hold on
p2=plot(CAsim(1:1200),hReleaseKm,'-k')
hold on
p3=plot(CAsim(1:1200),hReleaseSm,'-m')
xlabel('Crank Angle (^{o})')
ylabel('Heat Released (J)')
legend([p1(1) p2(1) p3(1)],'Mean Cycle','Knock Cycle','No-Knock Cycle','Location','northwest')
 
figure()
p1=plot(CAsim,Lossm,'-r')
hold on
p2=plot(CAsim,LossKm,'-k')
hold on
p3=plot(CAsim,LossSm,'-m')
xlabel('Crank Angle (^{o})')
ylabel('Heat Loss Rate (J/CA)')
legend([p1(1) p2(1) p3(1)],'Mean Cycle','Knock Cycle','No-Knock Cycle','Location','northeast')
 
figure()
p1=plot(CAsim,LossK,'-r','LineWidth',0.3) %T bulk S
hold on
p2=plot(CAsim,LossKm,'y','LineWidth',1.5)
hold on
xlabel('Crank Angle (^{o})')
ylabel('Heat Release Rates (J/(CA-m^{2}-K')
legend([p1(1) p2(1)],'All Cycles','Average Cycle')
xlim([0 30])
 
figure()
p1=plot(CAsim,LossS,'-r','LineWidth',0.3) %T bulk S
hold on
p2=plot(CAsim,LossSm,'y','LineWidth',1.5)
hold on
xlabel('Crank Angle (^{o})')
ylabel('Heat Release Rates (J/(CA-m^{2}-K')
legend([p1(1) p2(1)],'All Cycles','Average Cycle')
xlim([0 30])
 
figure()
p1=plot(CAsim,COVLossK,'-r','LineWidth',1.5) %T bulk S
hold on
p2=plot(CAsim,COVLossS,'-k','LineWidth',1.5)
hold on
xlabel('Crank Angle (^{o})')
ylabel('Coefficient of Variation')
legend([p1(1) p2(1)],'Knock','No-Knock','Location','northwest')
xlim([0 30])
 
figure()
p1=plot(CAsim(1:1200),MFBm,'-r')
hold on
p2=plot(CAsim(1:1200),MFBKm,'-k')
hold on
p3=plot(CAsim(1:1200),MFBSm,'-m')
xlabel('Crank Angle (^{o})')
ylabel('Mass Fraction Burned')
legend([p1(1) p2(1) p3(1)],'Mean Cycle','Knock Cycle','No-Knock Cycle','Location','northwest')
 
figure()
p1=plot(CAsim(1:1200),TotalHeatKm,'-r','LineWidth',1.5) %T bulk S
hold on
p2=plot(CAsim(1:1200),TotalHeatSm,'-g','LineWidth',1.5)
hold on
xlabel('Crank Angle (^{o})')
ylabel('Total Heat Released (J)')
legend([p1(1) p2(1)],'Knock','No-Knock','Location','northwest')
 
figure()
p1=plot(CAsim,TIKm,'-r','LineWidth',1.5) %T bulk S
hold on
p2=plot(CAsim,TISm,'-c','LineWidth',1.5)
hold on
xlabel('Crank Angle (^{o})')
ylabel('Isentropic Unburned Temperature (K)')
legend([p1(1) p2(1)],'Knock','No-Knock','Location','northwest')
 
 
c=input('Would you like to see all plots? \nif yes enter 1\nif no enter 2\n')
if c==2 
    close all
end
 
 
