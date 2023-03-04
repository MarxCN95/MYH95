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
%% Simulate the schedule
dt = [0.001*86400;1*86400;10*86400;100*86400;100*86400;200*86400;1*year;1*year;1*year;1*year;1*year;1*year;1*year;1*year];%300*86400;300*86400;300*86400;300*86400;300*86400;300*86400;300*86400];%单位为秒，86400为一天。
schedule = simpleSchedule(dt, 'W', W);
pcswitch = 1;%PR-EOS考虑毛管力的开关
pcapillary = 0 ;

%% Pick linear solver

linsolve = GMRES_ILUSolverAD('tolerance', 1e-2, 'maxIterations', 100);
linsolve.reduceToCell = false;
% Instead, we only keep the cell variables using the low-level interface.
linsolve.keepNumber = (ncomp+model.water)*G.cells.num;
% Right diagonal scaling helps with the condition number of the system due
% to mixed derivatives
linsolve.applyRightDiagonalScaling = true;
nls = NonLinearSolver('LinearSolver', linsolve);
[ws, states, rep] = simulateScheduleAD(state0, model, schedule);
%[ws, states, rep] = simulateScheduleAD(state0, model, schedule, 'nonlinearsolver', nls);
%[wsDiagonal, statesDiagonal, repDiagonal] = simulateScheduleAD(state0, modelDiagonalAD, schedule, 'nonlinearsolver', nls);
%%
getTime = @(report) [sum(cellfun(@(x) x.AssemblyTime, report.ControlstepReports{1}.StepReports{1}.NonlinearReport)), ... % Assembly
                     sum(cellfun(@(x) x.LinearSolver.SolverTime, report.ControlstepReports{1}.StepReports{1}.NonlinearReport(1:end-1))),... % Linear solver
                     sum(report.SimulationTime), ...
                     ];
time = getTime(rep);
%time_sparse = getTime(rep);
%time_Diagonal = getTime(repDiagonal);
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
