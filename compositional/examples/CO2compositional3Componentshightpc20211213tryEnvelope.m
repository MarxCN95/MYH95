%% -------------------------------------------含有毛管力----------------------------------------------------------
clear
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm linearsolvers
global pcswitch fluidint pcapillary;
global radius;
radius = 10;
includeWater = true;% Include aqueous phase
%includeWater = false;%true;% Don't Include aqueous phase


%[filename,pathname]=uigetfile('.data','请打开data文件');
%fn = fullfile(pathname,filename);%//修改，选择文件
%fn    = fullfile('C:','Users','MYH','Desktop','YP2composition','YP2E300v1','YP2E300_E300.DATA');
fn    = fullfile('C:','Users','MYH','Desktop','YP2composition','YP2E300v2','YP2E300V2_E300.DATA');
%fn    = fullfile('E:','Research','MRSTmodel','big grid3','1015_E300.DATA');
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck); 

G = initEclipseGrid(deck);
G = computeGeometry(G);
[Dx,Dy,Dz]=deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
%% Set up initial Fluid
fluidint = initDeckADIFluid(deck);
fluidint.rhoWS = 1000;
fluidint.rhoOS = 800;
fluidint.rhoGS = 100;
%fluidint.muW = @(p, varargin) constantViscosity(opt.mu(i), p, [1, 1,1]*centi*poise);%data里有pvtw就不需要设置了

eos = initDeckEOSModel(deck);
r = radius;
[Tcp, Pcp, deltaT, deltaP] = computeShift(eos.fluid, r);
%[Tcpi,Pcpi,deltaTi,deltaPi]= computedeltaShift(eos.fluid);
Shiftswitch = 0;
%Shiftswitch = true;
if Shiftswitch==1
    eos.fluid.Tcrit = Tcp;
    eos.fluid.Pcrit = Pcp;
end     
model = NaturalVariablesCompositionalModel(G, rock, fluidint, eos.fluid, 'water', includeWater);
%modelDiagonalAD = NaturalVariablesCompositionalModel(G, rock, fluidint, eos.fluid, 'water', includeWater, 'AutoDiffBackend', DiagonalAutoDiffBackend('modifyOperators', true));


%% Set up well pattern
pv = poreVolume(G, rock);
totTime = 9*year;
injrate = sum(pv)/totTime;
injrate = 0;
%Injection is pure gas
W = [];
 % Injector
    W = verticalWell(W, G, rock, 1, 1, [], 'comp_i', [0, 0, 1], ...
        'Type', 'rate', 'name', 'Inj', 'Val', injrate, 'sign', 1);

      for i = 1:numel(W)
         W(i).components = [1, 0, 0, 0, 0, 0,0];
     end
%% Set up reservior condition
ncomp = eos.fluid.getNumberOfComponents();
%6个组分
z0 = [0.012,0.387,0.194,0.119,0.235,0.041,0.012];
T = 0 + 273.15;
for j=1:13
p = 100*barsa;    
%figure;
for i=1:1000
%% Set up sat
if includeWater
    s0 = [0.35, 0.65, 0];
else
    s0 = [1, 0];
    for i = 1:numel(W)
        W(i).compi = W(i).compi(2:end);
    end
end
%% Set up initial state

pcswitch = 0;
state0 = initCompositionalState(G, p, T, s0, z0, eos);
%% Set up schedule
dt = [0.00001*86400];%0.2*86400;0.2*86400;0.2*86400;0.2*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400];%1*year;300*86400;300*86400;300*86400;300*86400;300*86400;300*86400;300*86400];%单位为秒，86400为一天。
schedule = simpleSchedule(dt, 'W', W);


%% start simulation
pcswitch = 1;%PR-EOS考虑毛管力的开关
pcapillary = 0 ;
[ws, states, rep] = simulateScheduleAD(state0, model, schedule);
L = states{1, 1}.L(36,1);
 value = (L-1)^2;
   if (L<1)&& value<0.1
        bubble(j,3)=L;
        bubble(j,2)=p;
        bubble(j,1)=T;
        save bubble;
        break
    end
%     if L>=1
%         p = (90*barsa+p)/2;
%     end
%     if L<1
%         p = (170*barsa+p)/2;
%     end
p = p-(0.1*barsa);
end
T = T+10;
end


pcswitch = 1;%PR-EOS考虑毛管力的开关
pcapillary = 500 ;
[wspc, statespc, reppc] = simulateScheduleAD(state0, model, schedule);

getTime = @(report) [sum(cellfun(@(x) x.AssemblyTime, report.ControlstepReports{1}.StepReports{1}.NonlinearReport)), ... % Assembly
                     sum(cellfun(@(x) x.LinearSolver.SolverTime, report.ControlstepReports{1}.StepReports{1}.NonlinearReport(1:end-1))),... % Linear solver
                     sum(report.SimulationTime), ...
                     ];
time = getTime(rep);

%  [schstep, qOs_D, qGs_D, qGs_T, qWs_D, qLs_D, fw, qOs_T, ...
%      qWs_T, qLs_T, bhp] = GeneratePRO(ws, schedule, 'PROD', 'Units', 'Field');
%% Plot well results
% plotWellSols(ws, dt);%绘制井曲线
% plotWellSols(wspc, dt);
% %% Plot the results in the interactive viewer
% v = [30, 60];
% figure('name','不考虑毛管力','Position',[100, 300, 1000, 400]); %clf;
% plotToolbar(G, states)
% view(v);
% axis tight
% colorbar('vert')
% 
% 
% v = [30, 60];
% figure('name','考虑毛管力','Position',[100, 300, 1000, 400]); %clf;
% plotToolbar(G, statespc)
% view(v);
% axis tight
% colorbar('vert')
% %% 液相组分含量
% figure;
% compliquidM = [ ] ;
% for step = 1:numel(statespc)
%     for i = 1:ncomp
%     compliquidM(step,i) = statespc{step,1}.x(36,i);
%     end
% end
% plot(dt,compliquidM);
% legend(eos.fluid.names);
% %% %绘制states.y，基质中各组分的气相百分比。
% figure;
% compvaporM = [ ] ;
% for step = 1:numel(statespc)
%     for i = 1:ncomp
%     compvaporM(step,i) = statespc{step,1}.y(36,i);
%     end
% end
% plot(dt,compvaporM);
% legend('气相', eos.fluid.names);