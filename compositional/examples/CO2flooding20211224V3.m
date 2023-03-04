%% -------------------------------------------含有毛管力----------------------------------------------------------
%%  仅临界性质偏移2和双重效应3
clear
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm linearsolvers
global pcswitch fluidint pcapillary;
global radius;
radius = 10;
includeWater = true;% Include aqueous phase
%includeWater = false;%true;% Don't Include aqueous phase
Shift = 1;  %critical properties shift

%[filename,pathname]=uigetfile('.data','请打开data文件');
%fn = fullfile(pathname,filename);%//修改，选择文件
%fn    = fullfile('C:','Users','MYH','Desktop','YP2composition','YP2E300v1','YP2E300_E300.DATA');
%fn    = fullfile('C:','Users','MYH','Desktop','YP2composition','YP2E300v2','YP2E300V2_E300.DATA');
%fn    = fullfile('C:','Users','MYH','Desktop','YP2composition','YP2E300v3fracturev3','YP2E300V3_E300.DATA');
%fn    = fullfile('E:','Research','MRSTmodel','big grid3','1015_E300.DATA');
fn    = fullfile('E:','YP2composition','YP2E300v3fracturev3','YP2E300V3_E300.DATA');
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck); 
% deck.GRID.PERMX = deck.GRID.PERMX*10;
% deck.GRID.PERMY = deck.GRID.PERMX*10;
% deck.GRID.PERMZ = deck.GRID.PERMX*10;
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

if Shift==1
    eos.fluid.Tcrit = Tcp;
    eos.fluid.Pcrit = Pcp;
end     
model = NaturalVariablesCompositionalModel(G, rock, fluidint, eos.fluid, 'water', includeWater);
%modelDiagonalAD = NaturalVariablesCompositionalModel(G, rock, fluidint, eos.fluid, 'water', includeWater, 'AutoDiffBackend', DiagonalAutoDiffBackend('modifyOperators', true));

%% Set up well pattern
pv = poreVolume(G, rock);
totTime = 1*year;
injrate = sum(pv)/totTime;
injrate =1*injrate;
%Injection is pure gas
%%
W = [];
 % Injector
    W = verticalWell(W, G, rock, 13, 1, [], 'comp_i', [0, 0, 1], ...
        'Type', 'rate', 'name', 'Inj1', 'Val', 0, 'sign', 1);
    W = verticalWell(W, G, rock, 25, 13, [], 'comp_i', [0, 0, 1], ...
        'Type', 'rate', 'name', 'Inj2', 'Val', 0, 'sign', 1);
    W = verticalWell(W, G, rock, 13, 25, [], 'comp_i', [0, 0, 1], ...
        'Type', 'rate', 'name', 'Inj3', 'Val', 0, 'sign', 1);
    W = verticalWell(W, G, rock, 1, 13, [], 'comp_i', [0, 0, 1], ...
        'Type', 'rate', 'name', 'Inj4', 'Val', 0, 'sign', 1);    
    W = verticalWell(W, G, rock, 13, 13, [],'comp_i', [0, 0.5, 0.5], ...
         'Name', 'Prod','Type', 'bhp', 'Val', 30*barsa, 'sign', -1);
      for i = 1:numel(W)
         W(i).components = [1, 0, 0, 0, 0, 0,0];
      end
        [Wprod, Wflooding] = deal(W);
        Wflooding(1).val   = injrate;
        Wflooding(2).val   = injrate;
        Wflooding(3).val   = injrate;
        Wflooding(4).val   = injrate;
%% 井坐标
% Wcell = [];
%  % Injector
%     Winjcell = verticalWell(Wcell, G, rock, 1, 1, [], 'comp_i', [0, 0, 1], ...
%         'Type', 'rate', 'name', 'Inj', 'Val', injrate, 'sign', 1);
% %      W = verticalWell(W, G, rock, 1, 1, [], 'comp_i', [0, 0, 1], ...
% %          'name', 'Inj', 'Type', 'bhp', 'Val', 200*barsa,'sign', 1);
% 
%     Wprodcell = verticalWell(Wcell, G, rock, Dx, Dy, [],'comp_i', [0, 0.5, 0.5], ...
%          'Name', 'Prod','Type', 'bhp', 'Val', 30*barsa, 'sign', -1);
%       for i = 1:numel(Wcell)
%          Wcell(i).components = [1, 0, 0, 0, 0, 0,0];
%      end
% %% 井设置
%   W = [];
%     Winj = addWell(W, G, rock, Winjcell.cells , 'comp_i', [0, 0, 1], ...
%         'Type', 'rate', 'name', 'Inj', 'Val', injrate, 'sign', 1);
%    Winj0 = addWell(W, G, rock, Winjcell.cells , 'comp_i', [0, 0, 1], ...
%         'Type', 'rate', 'name', 'Inj', 'Val', 0, 'sign', 1);
%  for i = 1:numel(Winj)
%     Winj(i).components = [1, 0, 0, 0, 0, 0,0];
%  end
%  for i = 1:numel(Winj0)
%     Winj0(i).components = [1, 0, 0, 0, 0, 0,0];
% end 
%     Wprod = addWell(W, G, rock, Wprodcell.cells,'comp_i', [0, 0.5, 0.5], ...
%           'Name', 'Prod','Type', 'bhp', 'Val', 30*barsa, 'sign', -1);
%  for i = 1:numel(W)
%     Wprod(i).components = [1, 0, 0, 0, 0, 0,0];
% end
%     
%% Set up reservior condition
ncomp = eos.fluid.getNumberOfComponents();
%6个组分
z0 = [0.012,0.387,0.194,0.119,0.235,0.041,0.012];
T = 60 + 273.15;
p = 150*barsa;
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
%% 断点继续算
% state0.pressure = states{34, 1}.pressure;
% state0.flux = states{34, 1}.flux;
% state0.s = states{34, 1}.s;
% state0.T = states{34, 1}.T;
% state0.components = states{34, 1}.components;
% state0.L = states{34, 1}.L;
% state0.K = states{34, 1}.K;
% state0.Z_V = states{34, 1}.Z_V;
% state0.x = states{34, 1}.x;
% state0.y = states{34, 1}.y;
% state0.Z_L= states{34, 1}.Z_L;
%% Set up schedule
dt = [0.001*86400;10*86400;30*86400;30*86400;30*86400;100*86400;100*86400;300*86400;300*86400;300*86400;300*86400;10*86400;30*86400;100*86400;100*86400];%100*86400;100*86400;100*86400;100*86400;100*86400;100*86400];%1*year;300*86400;300*86400;300*86400;300*86400;300*86400;300*86400;300*86400];%单位为秒，86400为一天。
%      1           2         3        4       5        6         7 1year      8%       9         10        11         12    13       14        15       
%schedule = simpleSchedule(dt, 'W', W);
dt = [];
dt = [0.001*86400;1*86400;10*86400;30*86400;100*86400];
for ti=6:9
    dt(ti) =1*year;
end
for ti=10:60
    dt(ti) =30*86400;
end
%dt = [0.001*86400;10*86400;30*86400;100*86400;1*year;1*year;1*year;1*year;1*year;1*year;1*year;1*year;1*year;1*year;1*year;1*year;1*year;1*year;1*year;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400;30*86400];%单位为秒，86400为一天。

    schedule = struct();
    schedule.control = [struct('W', Wprod);... % Production control 1
                        struct('W', Wflooding)]; % CO2 flooding control 2
 dtnumb = numel(dt);
 schedule.step.val = dt;
 schedule.step.control = ones(dtnumb,1);
 schedule.step.control(10:end,1) = 2;          
%% 但一时间步，为了查看网格是否正确
% dt = [];
% dt = [0.001*86400];
%     schedule = struct();
%     schedule.control = [struct('W', Wprod);... % Production control 1
%                         struct('W', Wflooding)]; % CO2 flooding control 2
%  dtnumb = numel(dt);
%  schedule.step.val = dt;
%  schedule.step.control = ones(dtnumb,1);
%  schedule.step.control(1:end,1) = 2;
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
caxis([3000000 15000000])
%%
v = [30, 60];
figure('name','不考虑毛管力','Position',[100, 300, 1000, 400]); %clf;
plotToolbar(G, rock)
view(v);
axis tight
colorbar('vert')
%%
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