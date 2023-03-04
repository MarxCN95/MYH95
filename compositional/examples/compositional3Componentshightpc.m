%% -------------------------------------------含有毛管力----------------------------------------------------------
clear
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm
global pcswitch fluidint;
%[filename,pathname]=uigetfile('.data','请打开data文件');
%fn = fullfile(pathname,filename);%//修改，选择文件
fn    = fullfile('C:','Users','myh','Desktop','0627','0623_E300.DATA');

deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
%%
fluidint = initDeckADIFluid(deck);

fluidint.rhoOS = 800;
fluidint.rhoGS = 10;

eos = initDeckEOSModel(deck);

model = NaturalVariablesCompositionalModel(G, rock, fluidint, eos.fluid, 'water', false);

schedule = convertDeckScheduleToMRST(model, deck);
[schedule.control.W.components] = deal([0.5, 0.5, 0, 0, 0, 0, 0]);%7个组分

%Injection is pure gas
%[schedule.control.W.compi] = deal([1, 0]);
%% Set up initial state
ncomp = eos.fluid.getNumberOfComponents();
%7个组分
z0 = [0.6, 0.3, 0.1, 0, 0, 0, 0];
T = 60 + 273.15;
p = 500*barsa;

pcswitch = 0;
state0 = initCompositionalState(G, p, T, [0.5, 0.5], z0, eos);
%% Simulate the schedule
pcswitch = 1;%PR-EOS考虑毛管力的开关
[ws, states, rep] = simulateScheduleAD(state0, model, schedule);
%%
[schstep, qOs_D, qGs_D, qGs_T, qWs_D, qLs_D, fw, qOs_T, ...
    qWs_T, qLs_T, bhp] = GeneratePRO(ws, schedule, 'PROD', 'Units', 'Field');
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
        title(eos.fluid.names{i})
        caxis([0, 1])
        colorbar('vert') %增加一个垂直色轴
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
figure(1); %clf;
plotToolbar(G, states)
view(v);
axis tight
colorbar('vert')
close all;
%% %绘制states.pressure，裂缝压力。
figure;
x=1:numel(states);
statespressure = [ ];
for step = 1:numel(states)
    statespressure(step,1) = states{step,1}.pressure(50,1);
end
plot(schstep,statespressure);
legend('fracturepressure');
%% %绘制wsgor，GOR。
figure;
x=1:numel(states);
wsgor = [ ] ;
for step = 1:numel(states)
    wsgor(step,1) = ws{step,1}.gor;
end
plot(schstep,wsgor);
legend('GOR');
%% %绘制states.rho，油相密度。
figure;
x=1:numel(states);
statesoilrho = [ ] ;
for step = 1:numel(states)
    statesoilrho(step,1) = states{step,1}.rho(50,1);
end
plot(schstep,statesoilrho);
legend('oilrho');
close all;
%%

%% --------------------------------------------------无毛管力------------------------------------------------------
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm
global pcswitch fluidint;
%[filename,pathname]=uigetfile('.data','请打开data文件');
%fn = fullfile(pathname,filename);%//修改，选择文件
fnnon    = fullfile('C:','Users','myh','Desktop','0623','0623_E300.DATA');

decknon = readEclipseDeck(fnnon);
decknon = convertDeckUnits(deck);

Gnon = initEclipseGrid(decknon);
Gnon = computeGeometry(Gnon);

rocknon  = initEclipseRock(decknon);
rocknon  = compressRock(rocknon, Gnon.cells.indexMap);
%%
fluidintnon = initDeckADIFluid(decknon);

fluidintnon.rhoOS = 800;
fluidintnon.rhoGS = 10;

eosnon = initDeckEOSModel(decknon);

modelnon = NaturalVariablesCompositionalModel(Gnon, rocknon, fluidintnon, eosnon.fluid, 'water', false);

schedulenon = convertDeckScheduleToMRST(modelnon, decknon);

[schedulenon.control.W.components] = deal([1, 0, 0, 0, 0, 0, 0]);%7个组分

%Injection is pure gas
%[schedule.control.W.compi] = deal([1, 0]);

%% Set up initial state
ncompnon = eosnon.fluid.getNumberOfComponents();
%7个组分
z0 = [0.6, 0.3, 0.1, 0, 0, 0, 0];
T = 60 + 273.15;
p = 500*barsa;
pcswitch = 0;
state0non = initCompositionalState(Gnon, p, T, [0.5, 0.5], z0, eosnon);
%% Simulate the schedule
pcswitch = 0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[wsnon, statesnon, repnon] = simulateScheduleAD(state0non, modelnon, schedulenon);
%%
[schstep, qOs_Dnon, qGs_Dnon, qGs_Tnon, qWs_Dnon, qLs_Dnon, fwnon, qOs_Tnon, ...
    qWs_Tnon, qLs_Tnon, bhpnon] = GeneratePRO(wsnon, schedulenon, 'PROD', 'Units', 'Field');
%% plot figure
figure(2); %clf;
plotToolbar(Gnon, statesnon)
view(v);
axis tight
colorbar('vert')
%% %绘制states.x，裂缝中各组分的液相百分比。
figure;
compliquid = [ ] ;
for step = 1:numel(states)
    for i = 1:ncomp
    compliquid(step,i) = states{step,1}.x(50,i);
    end
end
plot(schstep,compliquid);
legend('液相', eos.fluid.names);
%% %绘制states.x，基质中各组分的液相百分比。
figure;
compliquidM = [ ] ;
for step = 1:numel(states)
    for i = 1:ncomp
    compliquidM(step,i) = states{step,1}.x(52,i);
    end
end
plot(schstep,compliquidM);
legend(eos.fluid.names);

%% %绘制statesnon.x，裂缝中各组分的液相百分比。
figure;
compliquidnon = [ ] ;
for step = 1:numel(statesnon)
    for i = 1:ncomp
    compliquidnon(step,i) = statesnon{step,1}.x(50,i);
    end
end
plot(schstep,compliquidnon);
legend('液相', eosnon.fluid.names);
%% %绘制statesnon.x，基质中各组分的液相百分比。
figure;
compliquidnonM = [ ] ;
for step = 1:numel(statesnon)
    for i = 1:ncomp
    compliquidnonM(step,i) = statesnon{step,1}.x(52,i);
    end
end
plot(schstep,compliquidnonM);
legend('液相', eosnon.fluid.names);
%% %绘制states.y，裂缝中各组分的气相百分比。
figure;
compvapor = [ ] ;
for step = 1:numel(states)
    for i = 1:ncomp
    compvapor(step,i) = states{step,1}.y(50,i);
    end
end
plot(schstep,compvapor);
legend('气相', eos.fluid.names);
%% %绘制states.y，基质中各组分的气相百分比。
figure;
compvaporM = [ ] ;
for step = 1:numel(states)
    for i = 1:ncomp
    compvaporM(step,i) = states{step,1}.y(52,i);
    end
end
plot(schstep,compvaporM);
legend('气相', eos.fluid.names);
%% %绘制statesnon.y，裂缝中各组分的气相百分比。
figure;
compvapornon = [ ] ;
for step = 1:numel(statesnon)
    for i = 1:ncomp
    compvapornon(step,i) = statesnon{step,1}.y(50,i);
    end
end
plot(schstep,compvapornon);
legend('气相', eosnon.fluid.names);
%% %绘制statesnon.y，基质中各组分的气相百分比。
figure;
compvapornonM = [ ] ;
for step = 1:numel(statesnon)
    for i = 1:ncomp
    compvapornonM(step,i) = statesnon{step,1}.y(52,i);
    end
end
plot(schstep,compvapornonM);
legend('气相', eosnon.fluid.names);
%% %绘制states.components，裂缝中各组分的摩尔百分比。
figure;
x=1:numel(states);
statescomponents = [ ] ;
for step = 1:numel(states)
    for i = 1:ncomp
    statescomponents(step,i) = states{step,1}.components(50,i);
    end
end
plot(schstep,statescomponents);
legend(eos.fluid.names);
%% %绘制states.components，基质中各组分的摩尔百分比。
figure;
x=1:numel(states);
statescomponentsM = [ ] ;
for step = 1:numel(states)
    for i = 1:ncomp
    statescomponentsM(step,i) = states{step,1}.components(52,i);
    end
end
plot(schstep,statescomponentsM);
legend(eos.fluid.names);
%% %绘制statesnon.components，裂缝中各组分的摩尔百分比。
figure;
statescomponentsnon = [ ] ;
for step = 1:numel(statesnon)
    for i = 1:ncomp
    statescomponentsnon(step,i) = statesnon{step,1}.components(50,i);
    end
end
plot(schstep,statescomponentsnon);
legend(eosnon.fluid.names);
%% %绘制statesnon.components，基质中各组分的摩尔百分比。
figure;
statescomponentsnonM = [ ] ;
for step = 1:numel(statesnon)
    for i = 1:ncomp
    statescomponentsnonM(step,i) = statesnon{step,1}.components(52,i);
    end
end
plot(schstep,statescomponentsnonM);
legend(eosnon.fluid.names);
%% %绘制states.pressure，裂缝压力。
figure;
statespressurenon = [ ];
for step = 1:numel(statesnon)
    statespressurenon(step,1) = statesnon{step,1}.pressure(50,1);
end
plot(schstep,statespressurenon);
legend('fracturepressure');
%% %绘制wsgornon，GOR。
figure;
wsgornon = [ ] ;
for step = 1:numel(statesnon)
    wsgornon(step,1) = wsnon{step,1}.gor;
end
plot(schstep,wsgornon);
legend('GORnon');
%% %绘制statesnon.rho，油相密度。
figure;
statesoilrhonon = [ ] ;
for step = 1:numel(statesnon)
    statesoilrhonon(step,1) = statesnon{step,1}.rho(50,1);
end
plot(schstep,statesoilrhonon);
legend('oilrhonon');
close all;

%% 绘制compareqOs_T
figure;
compareqOs = [ ] ;
for step = 1:numel(statesnon)
    for i = 1:2
        if i == 1 ;
    compareqOs(step,i) = qOs_T(step,1);
        end
        if i == 2 ;
    compareqOs(step,i) = qOs_Tnon(step,1);
        end
    end
end
plot(schstep,compareqOs);
legend('qOs_T','qOs_Tnon');
%% 绘制compareqGs_T
figure;
compareqGs = [ ] ;
for step = 1:numel(statesnon)
    for i = 1:2
        if i == 1 
    compareqGs(step,i) = qGs_T(step,1);
        end
        if i == 2 
    compareqGs(step,i) = qGs_Tnon(step,1);
        end
    end
end
plot(schstep,compareqGs);
legend('qGs_T','qGs_Tnon');
%% 绘制compareqgor
figure;
comparegor = [ ] ;
for step = 1:numel(statesnon)
    for i = 1:2
        if i == 1 
    comparegor(step,i) = ws{step,1}.gor;
        end
        if i == 2 
    comparegor(step,i) = wsnon{step,1}.gor;
        end
    end
end
plot(schstep,comparegor);
legend('GOR','GORnon');
%% 绘制comparepressure
figure;
comparepressure = [ ] ;
for step = 1:11
    for i = 1:2
        if i == 1 
    comparepressure(step,i) = states{2,1}.pressure(step,1);
        end
        if i == 2 
    comparepressure(step,i) = statesnon{2,1}.pressure(step,1);
        end
    end
end
x = [1,2,3,4,5,6,7,8,9,10,11];
plot(x,comparepressure);
legend('pressure','pressurenon');
%% %绘制states.s，基质中的含油饱和度。
figure;
statesSO = [ ] ;
for step = 1:numel(states)
    statesSO(step,1) = states{step,1}.s(50,1);
end
plot(schstep,statesSO);
legend('SO');
%% %绘制statesnon.s，裂缝中的含油饱和度。
figure;
statesSOM = [ ] ;
for step = 1:numel(states)
    statesSOM(step,1) = states{step,1}.s(52,1);
end
plot(schstep,statesSOM);
legend('SO');
%% %绘制statesnon.s，裂缝中的含油饱和度。
figure;
statesSOnon = [ ] ;
for step = 1:numel(statesnon)
    statesSOnon(step,1) = statesnon{step,1}.s(50,1);
end
plot(schstep,statesSOnon);
legend('SO');
%% %绘制statesnon.s，基质中的含油饱和度。
figure;
statesSOnonM = [ ] ;
for step = 1:numel(statesnon)
    statesSOnonM(step,1) = statesnon{step,1}.s(52,1);
end
plot(schstep,statesSOnonM);
legend('SO');
%% %绘制states.Orho，基质中的油相密度。
figure;
statesOrho = [ ] ;
for step = 1:numel(states)
    statesOrho(step,1) = states{step,1}.rho(50,1);
end
plot(schstep,statesOrho);
legend('Orho');
%% %绘制statesnon.Orho，裂缝中的油相密度。
figure;
statesOrhoM = [ ] ;
for step = 1:numel(states)
    statesOrhoM(step,1) = states{step,1}.rho(52,1);
end
plot(schstep,statesOrhoM);
legend('Orho');
%% %绘制statesnon.Orho，裂缝中的油相密度。
figure;
statesOrhonon = [ ] ;
for step = 1:numel(statesnon)
    statesOrhonon(step,1) = statesnon{step,1}.rho(50,1);
end
plot(schstep,statesOrhonon);
legend('Orho');
%% %绘制statesnon.Orho，基质中的油相密度。
figure;
statesOrhononM = [ ] ;
for step = 1:numel(statesnon)
    statesOrhononM(step,1) = statesnon{step,1}.rho(52,1);
end
plot(schstep,statesOrhononM);
legend('Orho');
%% %绘制states.Grho，基质中的气相密度。
figure;
statesGrho = [ ] ;
for step = 1:numel(states)
    statesGrho(step,1) = states{step,1}.rho(50,2);
end
plot(schstep,statesGrho);
legend('Grho');
%% %绘制statesnon.Grho，裂缝中的气相密度。
figure;
statesGrhoM = [ ] ;
for step = 1:numel(states)
    statesGrhoM(step,1) = states{step,1}.rho(52,2);
end
plot(schstep,statesGrhoM);
legend('Grho');
%% %绘制statesnon.Grho，裂缝中的气相密度。
figure; 
statesGrhonon = [ ] ;
for step = 1:numel(statesnon)
    statesGrhonon(step,1) = statesnon{step,1}.rho(50,2);
end
plot(schstep,statesGrhonon);
legend('Grho');
%% %绘制statesnon.Grho，基质中的气相密度。
figure;
statesGrhononM = [ ] ;
for step = 1:numel(statesnon)
    statesGrhononM(step,1) = statesnon{step,1}.rho(52,2);
end
plot(schstep,statesGrhononM);
legend('Grho');

%%
state.huatu = 1;
figure(1); %clf;
plotToolbar(G, state)
view(v);
axis tight
colorbar('vert')
close all;