%% -------------------------------------------含有毛管力----------------------------------------------------------
clear
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm linearsolvers
global pcswitch fluidint pcapillary;
global radius;
radius = 20;
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
Shift = 0;
Shift = true;
if Shift==1
    eos.fluid.Tcrit = Tcp;
    eos.fluid.Pcrit = Pcp;
end     
model = NaturalVariablesCompositionalModel(G, rock, fluidint, eos.fluid, 'water', includeWater);
%modelDiagonalAD = NaturalVariablesCompositionalModel(G, rock, fluidint, eos.fluid, 'water', includeWater, 'AutoDiffBackend', DiagonalAutoDiffBackend('modifyOperators', true));


%% Set up well pattern
% pv = poreVolume(G, rock);
% totTime = 9*year;
% injrate = sum(pv)/totTime;
% 
% W = [];
% Winj = [];
% Wshut= [];
% % %% 定义井
% % % 定义生产井的位置
%   wellPoints = [5,10,4; 10,10,4; 15,10,4;20,10,4;25,10,4;30,10,4;35,10,4;40,10,4;45,10,4;50,10,4;55,10,4;60,10,4;65,10,4;70,10,4;75,10,4];
%  wellpath = makeSingleWellpath(wellPoints);
% % plotWellPath(wellpath);
% 
%   nperf = 75;%射孔总长度
%   perfstep = 5;%射孔间隔perfstep个网格
%   perfstartI = 5;%I方向射孔从perfstart网格开始
%   HI = (perfstartI :perfstep: nperf)';%水平井射孔的J网格
%   HJ = repmat(10, [nperf/perfstep, 1]);%水平井射孔的I网格
%   HK = repmat(5, [nperf/perfstep, 1]);%水平井射孔的K网格
%  wellPoints =[HI, HJ, HK];
%  wellpath = makeSingleWellpath(wellPoints);
% % plotWellPath(wellpath);
%   HcellInx = sub2ind(G.cartDims, HI, HJ, HK);% Convert IJK-indices to linear index (as used in G)
% 
%    W = addWell(W, G, rock, HcellInx, 'Type', 'bhp', ...
%       'Val', 60*barsa, 'Sign', -1, 'Comp_i', [0, 1, 0], 'Name', 'Producer1');
%   Winj = addWell(Winj, G, rock, HcellInx, 'Type', 'rate', ...
%     'Val', 0.0005, 'Sign', 1, 'Comp_i', [0, 0, 1], 'Name', 'Injector1');
% Wshut = addWell(Wshut, G, rock, HcellInx, 'Type', 'rate', ...
%     'Val', 0, 'Sign', 1, 'Comp_i', [0, 0, 1], 'Name', 'Injector1');
% 
% for i = 1:numel(W)
%     W(i).components = [1, 0, 0, 0, 0, 0,0];
% end
% for i = 1:numel(Winj)
%     Winj(i).components = [1, 0, 0, 0, 0, 0,0];
% end
% for i = 1:numel(Wshut)
%     Wshut(i).components = [1, 0, 0, 0, 0, 0,0];
% end

pv = poreVolume(G, rock);
totTime = 9*year;
injrate = sum(pv)/totTime;
%Injection is pure gas
W = [];
 % Injector
    W = verticalWell(W, G, rock, 1, 1, [], 'comp_i', [0, 0, 1], ...
        'Type', 'rate', 'name', 'Inj', 'Val', injrate, 'sign', 1);
%      W = verticalWell(W, G, rock, 1, 1, [], 'comp_i', [0, 0, 1], ...
%          'name', 'Inj', 'Type', 'bhp', 'Val', 200*barsa,'sign', 1);


    W = verticalWell(W, G, rock, Dx, Dy, [],'comp_i', [0, 0.5, 0.5], ...
         'Name', 'Prod','Type', 'bhp', 'Val', 10*barsa, 'sign', -1);
      for i = 1:numel(W)
         W(i).components = [1, 0, 0, 0, 0, 0,0];
     end
%% Set up reservior condition
ncomp = eos.fluid.getNumberOfComponents();
%6个组分
z0 = [0.012,0.387,0.194,0.119,0.235,0.041,0.012];
T = 60 + 273.15;
p = 170*barsa;
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
dt = [0.001*86400;1*86400;10*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400];%1*year;300*86400;300*86400;300*86400;300*86400;300*86400;300*86400;300*86400];%单位为秒，86400为一天。
schedule = simpleSchedule(dt, 'W', W);


% schedule = struct();
% 
% schedule.control(1) = struct('W', W); % productor control 1
% 
% schedule.control(2) = struct('W', Winj); % injector control 2
% 
% schedule.control(3) = struct('W', Wshut); % injector control 3
% 
% %schedule.control(4) = struct('W', W); % injector control 4
% 
% 
% dtnumb = numel(dt);
% schedule.step.val = dt;
% schedule.step.control = ones(dtnumb,1);
% schedule.step.control(7:9,1) = 2;
% schedule.step.control(10,1) = 3;
% %schedule.step.control(11:end,1) = 4;

%% start simulation
pcswitch = 0;%PR-EOS考虑毛管力的开关
pcapillary = 0 ;
[ws, states, rep] = simulateScheduleAD(state0, model, schedule);
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
plotWellSols(ws, dt);%绘制井曲线
plotWellSols(wspc, dt);
%% Plot the results in the interactive viewer
v = [30, 60];
figure('name','不考虑毛管力','Position',[100, 300, 1000, 400]); %clf;
plotToolbar(G, states)
view(v);
axis tight
colorbar('vert')


v = [30, 60];
figure('name','考虑毛管力','Position',[100, 300, 1000, 400]); %clf;
plotToolbar(G, statespc)
view(v);
axis tight
colorbar('vert')
%% 液相组分含量
figure;
compliquidM = [ ] ;
for step = 1:numel(statespc)
    for i = 1:ncomp
    compliquidM(step,i) = statespc{step,1}.x(36,i);
    end
end
plot(dt,compliquidM);
legend(eos.fluid.names);
%% %绘制states.y，基质中各组分的气相百分比。
figure;
compvaporM = [ ] ;
for step = 1:numel(statespc)
    for i = 1:ncomp
    compvaporM(step,i) = statespc{step,1}.y(36,i);
    end
end
plot(dt,compvaporM);
legend('气相', eos.fluid.names);