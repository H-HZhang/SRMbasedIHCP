function [ ReNewHeatFlux ] = Update3HeatFlux(kxls,Jacm,CalcTempdFT,qtrialFunc,nglcPara000,regParameter)
% qtrialFunc --- an a priori estimate of Q which can be set to zero if no a priori information is available.
% 2022/7/5
global nz2x NumRspTC RspTCNumLoc H123 Tsample Futuretime ;

% Iunit is the unit matrix for the regularization of heat flux during two
%% identity matrix
IdenMat =eye(nz2x); 
%% the sum of X'X for each Futuretime
SUM_XtX=zeros(nz2x,nz2x);
%% the sum of X'(Y-T) for each Futuretime
SUM_Xt_Y_T=zeros(nz2x,1);
for FtNum=1:Futuretime
    Xj=zeros(NumRspTC,nz2x);
    DetalT=zeros(NumRspTC,1);
    Xj =reshape(Jacm(:,FtNum),[NumRspTC nz2x]);
    [Y]=ReadTCdata(kxls,FtNum,Tsample,'interior');
    for ii=1:NumRspTC
        DetalT(ii,1)=Y(1,ii)-CalcTempdFT(RspTCNumLoc(ii,1),RspTCNumLoc(ii,2),FtNum);
    end
    SUM_XtX =SUM_XtX +Xj'*Xj;
    SUM_Xt_Y_T =SUM_Xt_Y_T +Xj'*DetalT;
end
%% the coefficients matrix of New Heat Flux
Argt=zeros(nz2x,nz2x);
% the coefficients matrix of Old Heat Flux
Blft=zeros(nz2x,1);
Argt=SUM_XtX+nglcPara000*IdenMat+regParameter*H123;
Blft=SUM_Xt_Y_T+(SUM_XtX+nglcPara000*IdenMat)*qtrialFunc;
% ReNewHeatFlux=Argt\Blft; %solve equ Argt*ANewHeatFlux=Blft
 ReNewHeatFlux = Argt\SUM_Xt_Y_T; %  a good choice
end

