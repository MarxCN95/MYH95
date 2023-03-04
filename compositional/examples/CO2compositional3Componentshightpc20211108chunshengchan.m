%% -------------------------------------------含有毛管力----------------------------------------------------------
clear
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm linearsolvers
global pcswitch fluidint pcapillary;

includeWater = true;% Include aqueous phase
%includeWater = false;%true;% Don't Include aqueous phase


%[filename,pathname]=uigetfile('.data','请打开data文件');
%fn = fullfile(pathname,filename);%//修改，选择文件
fn    = fullfile('C:','Users','MYH','Desktop','YP2composition','YP2E300v1','YP2E300_E300.DATA');
%fn    = fullfile('E:','Research','MRSTmodel','big grid3','1015_E300.DATA');
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck); 
%  deck.GRID.PERMX = deck.GRID.PERMX*100;
%  deck.GRID.PERMY = deck.GRID.PERMX*100;
%  deck.GRID.PERMZ = deck.GRID.PERMX*1000;
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

model = NaturalVariablesCompositionalModel(G, rock, fluidint, eos.fluid, 'water', includeWater);
modelDiagonalAD = NaturalVariablesCompositionalModel(G, rock, fluidint, eos.fluid, 'water', includeWater, 'AutoDiffBackend', DiagonalAutoDiffBackend('modifyOperators', true));
%schedule = convertDeckScheduleToMRST(model, deck);
%[schedule.control.W.components] = deal([1, 0, 0, 0, 0, 0]);%7个组分

%% Set up well pattern
pv = poreVolume(G, rock);
totTime = 9*year;
injrate = sum(pv)/totTime;
%Injection is pure gas
W = [];
% %% 定义井
% % 定义生产井的位置
  wellPoints = [5,10,4; 10,10,4; 15,10,4;20,10,4;25,10,4;30,10,4;35,10,4;40,10,4;45,10,4;50,10,4;55,10,4;60,10,4;65,10,4;70,10,4;75,10,4];
 wellpath = makeSingleWellpath(wellPoints);
 plotWellPath(wellpath);
% W = getWellFromPath([], G, rock, wellpath);
%  W = [];
  nperf = 75;%射孔总长度
  perfstep = 5;%射孔间隔perfstep个网格
  perfstartI = 5;%I方向射孔从perfstart网格开始
  HI = (perfstartI :perfstep: nperf)';%水平井射孔的J网格
  HJ = repmat(10, [nperf/perfstep, 1]);%水平井射孔的I网格
  HK = repmat(5, [nperf/perfstep, 1]);%水平井射孔的K网格
    wellPoints =[HI, HJ, HK];
 wellpath = makeSingleWellpath(wellPoints);
 plotWellPath(wellpath);
  HcellInx = sub2ind(G.cartDims, HI, HJ, HK);% Convert IJK-indices to linear index (as used in G)
% % 
% % % 
%  HIinj = (1 : nperf).' + 20;%水平井射孔的J网格
%  HJinj = repmat(5, [nperf, 1]);%水平井射孔的I网格
%  HKinj = repmat(10, [nperf, 1]);%水平井射孔的K网格
%  % % Convert IJK-indices to linear index (as used in G)
%  HcellInxinj = sub2ind(G.cartDims, HIinj, HJinj, HKinj);
   W = addWell(W, G, rock, HcellInx, 'Type', 'bhp', ...
      'Val', 10*barsa, 'Sign', -1, 'Comp_i', [0, 1, 0], 'Name', 'Producer');
%   Winj1 = addWell(W, G, rock, HcellInx, 'Type', 'rate', ...
%      'Val', 0.001, 'Sign', 1, 'Comp_i', [0, 0, 1], 'Name', 'Injector1');
%   W   = addWell(W, G, rock, HcellInxinj, 'Type', 'rate',...
%      'Val', 500/day, 'Sign',1, 'Comp_i', [1, 0], 'Name', 'Injector');
% %  plotGrid(G, 'FaceColor', 'g', 'FaceAlpha', .3, 'EdgeColor', 'w');
% % plotWell(G, W);
% % [nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));
% % cellinj = nx*ny*(nz-1)+1:nx:nx*ny*nz;
% % cellprod = nx:nx:nx*ny;
% % W   = addWell([], G, rock, flipud(cellinj), 'Type', 'rate',...
% %     'Val', 500/day, 'Sign',1, 'Comp_i', [1, 0], 'Name', 'Injector');
% % W   = addWell(W, G, rock, cellprod, 'Type', 'bhp', ...
% %     'Val', 100*barsa, 'Sign', -1, 'Comp_i', [0, 1], 'Name', 'Producer');
% plotGrid(G);
% view(30,50);
% plotWell(G,W);
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
%% Simulate the schedule
dt = [0.001*86400;1*86400;10*86400;100*86400;100*86400;200*86400;1*year;1*year;1*year;1*year;1*year;1*year;1*year];%300*86400;300*86400;300*86400;300*86400;300*86400;300*86400;300*86400];%单位为秒，86400为一天。
schedule = simpleSchedule(dt, 'W', W);
pcswitch = 1;%PR-EOS考虑毛管力的开关
pcapillary = 0 ;


%% start simulation

[ws, states, rep] = simulateScheduleAD(state0, model, schedule);

getTime = @(report) [sum(cellfun(@(x) x.AssemblyTime, report.ControlstepReports{1}.StepReports{1}.NonlinearReport)), ... % Assembly
                     sum(cellfun(@(x) x.LinearSolver.SolverTime, report.ControlstepReports{1}.StepReports{1}.NonlinearReport(1:end-1))),... % Linear solver
                     sum(report.SimulationTime), ...
                     ];
time = getTime(rep);

% [schstep, qOs_D, qGs_D, qGs_T, qWs_D, qLs_D, fw, qOs_T, ...
%     qWs_T, qLs_T, bhp] = GeneratePRO(ws, schedule, 'PROD', 'Units', 'Field');
%% Plot well results
plotWellSols(ws, dt);%绘制井曲线

%% Plot the results in the interactive viewer
v = [30, 60];
figure('name','不考虑毛管力','Position',[100, 400, 1300, 500]); %clf;
plotToolbar(G, states)
view(v);
axis tight
colorbar('vert')



% lf = get(0, 'DefaultFigurePosition');
% h = figure('Position', lf + [0, -200, 350, 200]);
% nm = ceil(ncomp/2);
% v = [-30, 60];
% for step = 1:numel(states)
%     figure(h); clf
%     state = states{step};
%     for i = 1:ncomp
%         subplot(nm, 3, i);
%         plotCellData(G, state.components(:, i), 'EdgeColor', 'none');
%         view(v);
%         title(eos.fluid.names{i})
%         caxis([0, 1])
%         colorbar('vert') %增加一个垂直色轴
%     end
%     subplot(nm, 3, ncomp + 1);
%     plotCellData(G, state.pressure, 'EdgeColor', 'none');
%     view(v);
%     title('Pressure')
%     
%     subplot(nm, 3, ncomp + 2);
%     plotCellData(G, state.s(:, 1), 'EdgeColor', 'none');
%     view(v);
%     title('sW')
%     
%     subplot(nm, 3, ncomp + 3);
%     plotCellData(G, state.s(:, 2), 'EdgeColor', 'none');
%     view(v);
%     title('sO')
%     drawnow
% end
