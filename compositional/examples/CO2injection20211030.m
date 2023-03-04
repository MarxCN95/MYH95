%% Example demonstrating a three dimensional, six component problem
% We set up a simple grid and rock structure. The numbers can be adjusted
% to get a bigger/smaller problem. The default is a small problem with
% 20x20x2 grid blocks.
clear all
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm
global pcswitch pcapillary

%fn = fullfile(pathname,filename);%//修改，选择文件
 fn    = fullfile('C:','Users','MYH','Desktop','1Dmodel','1DMODEL_E300.DATA');
% fn    = fullfile('E:','Research','MRSTmodel','big grid3','1015_E300.DATA');
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck); 

G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

%% Set up quarter five spot well pattern
% We place vertical wells in opposing corners of the reservoir. The
% injector is rate controlled and the producer is bottom hole pressure
% controlled.
W = [];
% Injector
% W = verticalWell(W, G, rock, 1, 1, [], 'comp_i', [1, 0], 'name', 'Inj',...
%     'Val', 0.0015, 'sign', 1, 'type', 'rate');
W = verticalWell(W, G, rock, 1, 1, [], 'comp_i', [1, 0], 'name', 'Inj',...
    'Val', 0.00005, 'sign', 1, 'type', 'rate');
% Producer
W = verticalWell(W, G, rock, 20, 1, [], ...
    'comp_i', [0.5, 0.5], 'Name', 'Prod', 'Val', 5000000, 'sign', -1, 'Type', 'bhp');
% %% 定义井
% % 定义生产井的位置
% % wellPoints = [328.55, 246.39; 412.19, 303.49; 547.63,395.39; 639.58,460.61; 1158.9,781.39; 1275.39, 865.53; 1394.86,937.67; 1493.71,1010.82; 1563.73,1054.44; 1608.14,1086.22; 1700.09,1128.99; 2141.58,1382.57];
% % wellpath = makeSingleWellpath(wellPoints);
% % plotWellPath(wellpath);
% % W = getWellFromPath([], G, rock, wellpath);
%  W = [];
%  nperf = 40;%射孔长度
%  HI = (1 : nperf).' + 20;%水平井射孔的J网格
%  HJ = repmat(1, [nperf, 1]);%水平井射孔的I网格
%  HK = repmat(5, [nperf, 1]);%水平井射孔的K网格
%  HcellInx = sub2ind(G.cartDims, HI, HJ, HK);
% % % Convert IJK-indices to linear index (as used in G)
% % % 
%  HIinj = (1 : nperf).' + 20;%水平井射孔的J网格
%  HJinj = repmat(5, [nperf, 1]);%水平井射孔的I网格
%  HKinj = repmat(10, [nperf, 1]);%水平井射孔的K网格
%  % % Convert IJK-indices to linear index (as used in G)
%  HcellInxinj = sub2ind(G.cartDims, HIinj, HJinj, HKinj);
%   W = addWell(W, G, rock, HcellInx, 'Type', 'bhp', ...
%      'Val', 100*barsa, 'Sign', -1, 'Comp_i', [0, 1], 'Name', 'Producer');
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
%%
fluidint = initDeckADIFluid(deck);

fluidint.rhoOS = 800;
fluidint.rhoGS = 10;

eos = initDeckEOSModel(deck);

model = NaturalVariablesCompositionalModel(G, rock, fluidint, eos.fluid, 'water', false);
%% Set up initial state
ncomp = eos.fluid.getNumberOfComponents();
%6个组分
z0 = [0.1, 0.3, 0.2, 0.2, 0.1, 0.1];
T = 60 + 273.15;
p = 120*barsa;

pcswitch = 0;
state0 = initCompositionalState(G, p, T, [0.5, 0.5], z0, eos);
 for i = 2:numel(W)
     W(i).components = [1,0,0,0,0,0];
 end
%% Set up schedule and simulate the problem
% We simulate two years of production with a geometric rampup in the
% timesteps.
%time = 10*year;
%n = 45;
%dt = rampupTimesteps(time, 100*day, 5);
dt = [0.001*86400;1*86400;10*86400;10*86400];%;100*86400;300*86400;100*86400;100*86400;100*86400];%单位为秒，86400为一天。
%schedule = simpleSchedule(dt, 'W', W);
%% Run 比较
pcswitch = 0;%%毛管力开关
[ws, states, rep] = simulateScheduleAD(state0, model, schedule);
pcswitch = 1;%%毛管力开关
pcapillary = 0.5*1000000 ;
[wscap, statescap, repcap] = simulateScheduleAD(state0, model, schedule);
%% Plot all the results
% lf = get(0, 'DefaultFigurePosition');
% h = figure('Position', lf + [0, -200, 350, 200]);
% nm = ceil(ncomp/2);
 v = [-30, 60];
% for step = 1:numel(states)
%     figure(h); clf
%     state = states{step};
%     for i = 1:ncomp
%         subplot(nm, 3, i);
%         plotCellData(G, state.components(:, i), 'EdgeColor', 'none');
%         view(v);
%         title(fluid.names{i})
%         caxis([0, 1])
%     end
%     subplot(nm, 3, ncomp + 1);
%     plotCellData(G, state.pressure, 'EdgeColor', 'none');
%     view(v);
%     title('Pressure')
%     
%     subplot(nm, 3, ncomp + 2);
%     plotCellData(G, state.s(:, 1), 'EdgeColor', 'none');
%     view(v);
%     title('sO')
%     
%     subplot(nm, 3, ncomp + 3);
%     plotCellData(G, state.s(:, 2), 'EdgeColor', 'none');
%     view(v);
%     title('sG')
%     drawnow
% end
%% Plot the results in the interactive viewer
set(0,'defaultfigurecolor','w') 
figure(1); clf;
plotToolbar(G, statespc)
view(v);
% axis tight
axis equal
caxis([1e6,1.2e7])
colorbar
hold on

caxis manual


 plotWellSols({ws,wscap}, cumsum(schedule.step.val), ...
     'datasetnames', {'NO', 'YES'}, 'linestyles', {'-', '-.'});
 box on
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
