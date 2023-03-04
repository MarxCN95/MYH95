%% Example demonstrating a three dimensional, six component problem
% We set up a simple grid and rock structure. The numbers can be adjusted
% to get a bigger/smaller problem. The default is a small problem with
% 20x20x2 grid blocks.
mrstModule add ad-core ad-props mrst-gui compositional
% Dimensions
nx = 20;
ny = nx;
nz = 2;
% Name of problem and pressure range
casename = 'lumped_1';
minP = 100*barsa;%//转换单位100bar=10MPa=10000000Pa,转换成pa
maxP = 200*barsa;%//转换单位
% Set up grid and rock
dims = [nx, ny, nz];
pdims = [1000, 1000, 1];
G = cartGrid(dims, pdims);%//创建一个20*20*2的网格，网格总的长宽高100*1000*1
G = computeGeometry(G);%//向网格中添加几何信息
rock = makeRock(G, 50*milli*darcy, 0.25);%//设置岩石渗透率50mD，孔隙度0.25
%% Set up quarter five spot well pattern
% We place vertical wells in opposing corners of the reservoir. The
% injector is rate controlled and the producer is bottom hole pressure
% controlled.%//直井布置在储层的角落。注入井由速率控制，生产井由井底压力控制。
W = [];
% Injector
W = verticalWell(W, G, rock, 1, 1, [], 'comp_i', [1, 0], 'name', 'Inj',...
    'Val', 0.0015, 'sign', 1, 'type', 'rate');
% Producer
W = verticalWell(W, G, rock, nx, ny, [], ...
    'comp_i', [0.5, 0.5], 'Name', 'Prod', 'Val', minP, 'sign', -1, 'Type', 'bhp');
%% Set up model and initial state
% We set up a problem with quadratic relative permeabilities. The fluid
% model is retrieved from "High Order Upwind Schemes for Two-Phase,
% Multicomponent Flow" (SPE 79691) by B. T. Mallison et al.
%
% The model consists of six components. Several of the components are not
% distinct molecules, but rather lumped groups of hydrocarbons with similar
% molecular weight. The reservoir contains all these components initially,
% which is then displaced by the injection of pure CO2.
%//建立模型和初始状态。我们建立了一个二次相对渗透率问题。流体模型由B. T. Mallison等人的“两相多组分流的高阶迎风格式”(SPE 79691)反演。
%//该模型由六个组分组成。有几个组分不是纯分子，而是分子量相似的碳氢化合物的集中基团。储层最初包含所有这些成分，然后被注入纯二氧化碳所取代。
nkr = 2;
[fluid, info] = getCompositionalFluidCase(casename);
flowfluid = initSimpleADIFluid('n', [nkr, nkr, nkr], 'rho', [1000, 800, 10]);

gravity reset on
model = NaturalVariablesCompositionalModel(G, rock, flowfluid, fluid, 'water', false);

ncomp = fluid.getNumberOfComponents();
state0 = initCompositionalState(G, minP + (maxP - minP)/2, info.temp, [0.5, 0.5], info.initial, model.EOSModel);

for i = 1:numel(W)
    W(i).components = info.injection;
end
%% Set up schedule and simulate the problem
% We simulate two years of production with a geometric rampup in the
% timesteps.
time = 3*day;
n = 45;
dt = rampupTimesteps(time, 1*day, 1);
schedule = simpleSchedule(dt, 'W', W);
pcswitch = 1;
[ws, states, rep] = simulateScheduleAD(state0, model, schedule);
%% %//
%sgof = cellfun(convert, slgof, 'UniformOutput', false);
%cfun  = @(f) cellfun(f, sgof, 'UniformOutput', false);
%Tpcog = extendTab( cfun(@(x) x(:, [1, 4])) );
% f.pcOG = @(sg, varargin) ireg(Tpcog, sg                   , varargin{:});
            
%            pc=cellfun(@(x)f.pcOG(x.s(:,2)),states,'UniformOutput',false);
%% Plot all the results
lf = get(0, 'DefaultFigurePosition');
h = figure('Position', lf + [0, -200, 350, 200]);
nm = ceil(ncomp/2);
v = [-30, 60];
for step = 1:numel(states)
    figure(h); clf
    state = states{step};
    for i = 1:ncomp
        subplot(nm, 3, i);
        plotCellData(G, state.components(:, i), 'EdgeColor', 'none');
        view(v);
        title(fluid.names{i})
        caxis([0, 1])
    end
    subplot(nm, 3, ncomp + 1);
    plotCellData(G, state.pressure, 'EdgeColor', 'none');
    view(v);
    title('Pressure')
    
    subplot(nm, 3, ncomp + 2);
    plotCellData(G, state.s(:, 1), 'EdgeColor', 'none');
    view(v);
    title('sO')
    
    subplot(nm, 3, ncomp + 3);
    plotCellData(G, state.s(:, 2), 'EdgeColor', 'none');
    view(v);
    title('sG')
    drawnow
end
%% Plot the results in the interactive viewer
figure(1); clf;
plotToolbar(G, states)
view(v);
axis tight
[schstep, qOs_D, qGs_D, qGs_T, qWs_D, qLs_D, fw, qOs_T, ...
    qWs_T, qLs_T, bhp] = GeneratePRO(ws, schedule, 'Prod', 'Units', 'Field');
%% Copyright notice
% <html>
% <p><font size="-1">
% Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
