% This is a 2D Heat Conduction Program for calculation of heat flux q2
% boundary condition. Firstly, a  Guess Heat Flux will apply on the hot
% side of computation area, where is correspond to the Mold Hot Surface.
%                       q1
%        j=nz    ----------------
%                | hot  .       | cold
%                |      .       |
%                |      .       |
%          q2?   | -3-  .       | Td
%                |      .       |
%                |      .       |
%                | face .       | face
%                |      .       |
%        j=1     ----------------
%               i=1     q3      i=nx
% 2022/7/3
% The Inverse Problem can be represented as an optimization problem
% the objective function
%       s = sum_Fts_(Y -T)^2 + beta*||h_dq||^2
% with partial differential equation (PDE) constraints
%       ∂T/∂t=ΔT with  T=T0, -∂T/∂n = q and T(Γ_4,t)=f(t)
% q - estimated one at time interval [tj,tj+tr], q - cal., one at latest time step tj-1
% for the time interval [tj,tj+tr], the q* is a time-independent and spatial-dependent varible.
% α - the regularization parameter (>0)
% h_d - dth-order derivative operator
% TcNumLoc_nznx - fetch the numeric loction for both cold side TC and RspTC.
% MatPropAreaDis - Material Property , calculated domain size and its discrete
% LoadTCdata - read tc data and Fs
% SolutionSetup - CalculateTime,nt,dt
% TcNumLoc_nznx - TC numeric location
% ReadTCdata - read Coldside/Rsptc temperature
%% clear workship
clear;close all;clc;
% declare the global variables
global H L kappa TCC heatcapacity density VolumetricHeat TC_Diameter_t Lref;
global Fs Tsample TCsample H123 NumColdTC NumRspTC RspTCNumLoc ColTCNumLoc res_mLoc cold_mLoc;
global dt CalculateTime nste Futuretime nx nz nz2x dx dz x2d z2d z_nz msh_num msh_x msh_z msh_bc;
% the global variables must be deleted before they could be used.
Tsample =[];NumColdTC=[]; NumRspTC=[]; RspTCNumLoc=[]; ColTCNumLoc=[];
res_mLoc =[]; cold_mLoc=[];x2d=[]; z2d=[];z_nz=[];
msh_num =[];msh_x=[];msh_z=[];msh_bc=[];
%% load data
disp('Attention data in LoadTCdata.m');
ExpNum=1;
dataName = 'CCmsExptData';
TC_Diameter_t = 0.0005; % 0.5 mm
[msh,Fst,nste,TCsample] = LoadEpxTcData( ExpNum,strcat(dataName,'.mat') );
TCsampleK = TCsample;
TCsampleK(:,2:end) = TCsample(:,2:end) +273.15;
%% Material Property , calculated domain size and its discrete(nz*nx ponits)
MolTyp=1;
[Ht,Lt,nx,nz,kappat,TCCt,heatcapacityt,densityt,VolumetricHeatt]=MatPropAreaDis(MolTyp,msh);
disp('Attention data in MatPropAreaDis.m');
%% Ref., Value and Normalization
Lref = Ht;qref = 1e6;Tmax_min= qref*Lref/TCCt;Tref =0;
H = Ht/Lref;
L = Lt/Lref;
kappa =1;
TCC=1;
heatcapacity=1;
density=1;
VolumetricHeat =heatcapacity*density;
% solve setup
Fs = 1/(1/Fst*kappat/Lref/Lref);
% Tsample =TCsample(:,1:17);
Tsample(:,1)= TCsample(:,1)*kappat/Lref/Lref;
Tsample(:,2:17)= (TCsample(:,2:17) -Tref)/Tmax_min;
%Std of Tc is 0.75%Tmax
normed_sigma2 = (0.0075*max(max(TCsample(:,2:17))) -Tref)/Tmax_min;

% TrueQ123 = TrueQ123t/qref;
disp('Attention data in SolutionSetup.m');
[CalculateTimet,dtt] = SolutionSetup( Fst,nste );
CalculateTime = CalculateTimet*kappat/Lref/Lref;
dt =dtt*kappat/Lref/Lref;
%meshing
nz2x=nz + nx+nx;
dx = L/(nx-1); % Spacing of grid in x-direction
dz = H/(nz-1); % Spacing of grid in z-direction
[x2d,z2d] = meshgrid(0:dx:L, 0:dz:H); % create even grid
z_nz=z2d(:,1);
%% mesh data
num = 1;
for i=1:nx
    for j=1:nz
        msh_num(i,j) = num;    %give each node a number
        msh_x(i,j) = (i-1)*dx;
        msh_z(i,j) = (j-1)*dz;
        num = num+1;
    end
end
% is point at boundary? 0- No, others- Yes.
msh_bc= zeros(size(msh_num));
msh_bc(1:nx,1) =1;
msh_bc(1:nx,nz) =1;
msh_bc(1,1:nz) =1;
msh_bc(nx,1:nz) =1;
%% TC numeric location
% *_mLoc = [x/mm, z/mm, the column of data in excel]
x_resT = 0.003;
res_mLoc = [
    x_resT 0.000 2;
    x_resT 0.003 4;
    x_resT 0.006 6;
    x_resT 0.009 8;
    x_resT 0.012 10;
    x_resT 0.015 12;
    x_resT 0.018 14;
    x_resT Ht    16;
    ];
cold_mLoc = [
    Lt 0.000 3;
    Lt 0.003 5;
    Lt 0.006 7;
    Lt 0.009 9;
    Lt 0.012 11;
    Lt 0.015 13;
    Lt 0.018 15;
    Lt Ht    17;
    ];
% *TCNumLoc = (numIdx-Z,numIdx-X,msh_num,the column of data in excel to be readed,dimensionless x,dimensionless z, at boundary?)
[RspTCNumLoc,x_Rspc,z_RspTc] =TCnumLoc(nx,nz,H,L,Lref,res_mLoc, 'interior');
[NumRspTC]=size(RspTCNumLoc,1);RspTCNumLocY=RspTCNumLoc(:,1);
[ColTCNumLoc,x_ColTc,z_ColTc] =TCnumLoc(nx,nz,H,L,Lref,cold_mLoc, 'coldside' );
[NumColdTC]=size(ColTCNumLoc,1);ColTCNumLocY=ColTCNumLoc(:,1);
Cstar = 709;
Cstop = 997;
%%
maxIteration = 1;
% total of the time steps
[tolTimeStep,~]=size(Tsample);
% Number of Future Time Steps
numFts= 6; %[1,2,4,6,8,10,12];
Pare0M = [0];           %, 0
% regularization parameter
regPara = 7.32e-5 ;% [1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,1e-1,1,5,10,30,40,75,100];
order_spatial_reg= {'zeros-spatial-reg', '1st-spatial-reg', '2nd-spatial-reg'};
currentFolder = pwd;
xLogYsbtT1 =zeros(length(regPara),length(numFts),length(order_spatial_reg))+NaN;
yLogLq1    =zeros(length(regPara),length(numFts),length(order_spatial_reg))+NaN;
photosfile_name =  strcat(currentFolder,'\','Pictures');
eval(['mkdir',32, photosfile_name]);

for id_order_spat_reg = 1  % 1:length(order_spatial_reg)
    spatial_reg = cell2mat(order_spatial_reg(id_order_spat_reg));
    [H_b1] = spaRegMatrix(nx,spatial_reg);
    [H_b2] = spaRegMatrix(nz,spatial_reg);
    [H_b3] = spaRegMatrix(nx,spatial_reg);
    % regularization matrix, H = h'*h
    H123 = blkdiag(H_b1,H_b2,H_b3);
    
    for idx_numFts = 1:length(numFts)
        Futuretime= numFts(idx_numFts);
        %  Jacobian matrix
        Jacm = JacobianMatrix3(0,nz2x,Futuretime);
        tol = norm( zeros(NumRspTC*(tolTimeStep-Futuretime),1) +  normed_sigma2);
        for i_Pare0 =1:length(Pare0M)
            
            Para0 = Pare0M(i_Pare0); % to be neglected parameter
            
            for idx_regPara=1:2:length(regPara)
                regPi = regPara(idx_regPara);
                %% make file
                str1 = strcat(spatial_reg,',Nqest msh,',num2str(msh),', nFts,',num2str(Futuretime),', beta,',num2str(regPi));
                file_name =  strcat(dataName,'-nFts',num2str(Futuretime),...
                    '-msh',num2str(msh),'-',spatial_reg,'-beta-',num2str(regPi));
                
                if exist(file_name,'dir') && exist(strcat(currentFolder,'\',file_name,'\',file_name,'.mat'))
                    load(strcat(currentFolder,'\',file_name,'\',file_name,'.mat'));
                else
                    eval(['mkdir',32,file_name]);
                    [~,Train_file_values] = fileattrib(file_name);
                    
                    % set initial temp at the beginning of calculation
                    [MolIntTemp] =InitialCondtion(1 ,Lref);
                    MolIntTemp = MolIntTemp*0 + mean(mean(MolIntTemp));
                    
                    MHFheatflux=zeros(tolTimeStep-Futuretime,nz2x);
                    AllTemp=zeros(nz,nx,tolTimeStep-Futuretime);
                    DiffErr=10;
                    RecRmsTstp =zeros(tolTimeStep-Futuretime-1,1);
                    RecRmsIter=sparse(tolTimeStep-Futuretime,maxIteration);
                    % start the timer
                    Lapsttime= 0;
                    tic;
                    for kxls=1:tolTimeStep-Futuretime       % time squence kxls
                        CrtTime= Tsample(kxls,1);  % current time
                        Q123 = zeros(nz2x,1);% nz2x*1 vector
                        
                        itern=1;
                        maxRms=0.005;
                        
                        %% the iterations
                        while DiffErr>=maxRms
                            Q3 = Q123(1:nx);
                            Q2 = Q123(nx+1:nx+nz);
                            Q1 = Q123(nx+nz+1:nz2x);
                            [CalTdFt,Timenow]=CNSolFourEquQ3(CrtTime,kxls,Futuretime,MolIntTemp,Q1,Q2,Q3);%CN
                            % Y-T
                            [ Y_T_nFts,Ytem,~] = resCalMea(kxls,Futuretime,Tsample,CalTdFt);
                            DiffErr= norm(Y_T_nFts)/norm(Ytem);
                            
                            if itern <= maxIteration
                                itern=itern+1;
                                Q123 = Update3HeatFlux(kxls,Jacm,CalTdFt,Q123,Para0,regPi);
                            else
                                break
                            end % end if
                            
                        end % end while
                        MHFheatflux(kxls,1:nz2x) = Q123(1:nz2x,1);
                        MolIntTemp=CalTdFt(:,:,1);
                        AllTemp(:,:,kxls)=CalTdFt(:,:,1);
                        DiffErr=10;
                    end % end kxls time squence
                    %                     disp(strcat(str ,' | Proceding:',num2str(kxls),'/',num2str(sm),' | maxit:',num2str(itern)));
                    % Train_file_status=0;
                    Lapsttime= toc;
                    
                end % exist(file_name,'dir')
                [YsbtT, Yture, Tcalu ] = resCalMea(1,tolTimeStep-Futuretime,Tsample,AllTemp);
                
                MHFQ2 = MHFheatflux(1:end,nx+1:nz+nx);
                if numel(MHFQ2(isnan(MHFQ2)))>1
                    disp(strcat(str1,', NaN in Q2'))
                else
                    normResYT = norm(YsbtT);
                    % L = { log(||Y -T||), log (||Lq||) } = {xlog_Y2T, ylog_Lq }
                    xLogYsbtT1(idx_regPara, idx_numFts, id_order_spat_reg) = log10( normResYT);
                    squLq=0;
                    for idxLq = 1:tolTimeStep-Futuretime
                        squLq = squLq + MHFheatflux(idxLq,:)*H123*MHFheatflux(idxLq,:)';
                    end
                    yLogLq1( idx_regPara, idx_numFts, id_order_spat_reg) = log10(sqrt(squLq));
                    
                    figure(1);clf;
                    plot(TCsample(1:tolTimeStep-Futuretime,1),MHFQ2);
                    
                    title({str1;strcat('tol=', num2str(tol), ' ||Y-T||=',num2str(normResYT),...
                        ' ||Lq||=',num2str(sqrt(squLq)),' regP*||Lq||^2=',num2str(regPi*squLq) ); strcat('Re_T=',num2str(100*normResYT/norm(Yture)),...
                        '% Re_e_x_p=', num2str(100*tol/norm(Yture)) ,'%')});
                    grid on;
                    %                     print(gcf,'-dpng',strcat(Train_file_values.Name,'\cost =',num2str(N2cost),'a','.png'));
                    print(gcf,'-dpng',strcat(photosfile_name,'\',file_name, '-norm(Y-T)=',num2str(normResYT),'.png'));
                    %                     OnlinePost(N2cost,MHFQ2, AllTemp*Tmax_min+Tref,Train_file_values.Name,709,997,Fst);
                    disp(strcat(str1,', tol,', num2str(tol), ', ||Y-T||,',num2str(normResYT),...
                        ', ||Lq||,',num2str(sqrt(squLq)), ', Re_T,',num2str(100*normResYT/norm(Yture)),...
                        '%, Re_exp,', num2str(100*tol/norm(Yture)) ,'%,'  ,file_name   ));
                end
                
                save(strcat(Train_file_values.Name,'\',file_name,'.mat'));
                
            end % end beta loop
        end % end alpha loop
    end % end r loop
end % end order loop
save(strcat('Res-',datestr(now,1), '.mat'  ));
%% data analysis
% select a result
file_name =  strcat(dataName,'-nFts',num2str(6),...
    '-msh',num2str(msh),'-',spatial_reg,'-beta-',num2str(regPi));
load(strcat(currentFolder,'\',file_name,'\',file_name,'.mat'));
% close all
[idx_Fig] = DrawFigures(normResYT,MHFQ2, AllTemp*Tmax_min+Tref,Train_file_values.Name,709,997,Fst);
[idx_Fig] = TempHeatFluxContour(AllTemp*Tmax_min+Tref+273.15,TCCt,340,390,...
    709,36,997,dx*Lref,dz*Lref,nz,nx,x2d*Lref,z2d*Lref,TCsample,idx_Fig);
print(gcf,'-dpng',strcat(photosfile_name,'\','Fig.',num2str(idx_Fig-1),'.png'));
%% Wavelet Transform for mold heat flux
figure(idx_Fig);
idx_Fig = idx_Fig +1;
cwt(MHFQ2(709:997,31),Fst);
figure(idx_Fig);
xlim([TCsample(709,1)  TCsample(997,1)]);
cwt(MHFQ2(:,31),seconds(1/Fst))
%  y=Acos(ωt+φ)+h, and cycle T=2π/ω=1/f
wcoherence(3*cos(2*pi*10/6*TCsample(1:tolTimeStep-Futuretime ,1)+ 60*pi/90),MHFQ2(:,31),Fst);
% xlim([11  16.5])
print(gcf,'-dpng',strcat(photosfile_name,'\','Fig.',num2str(idx_Fig-1),'.png'));
