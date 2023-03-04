%% -------------------------------------------����ë����----------------------------------------------------------
clear
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm
global pcswitch fluidint;
%[filename,pathname]=uigetfile('.data','���data�ļ�');
%fn = fullfile(pathname,filename);%//�޸ģ�ѡ���ļ�
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
[schedule.control.W.components] = deal([0.5, 0.5, 0, 0, 0, 0, 0]);%7�����

%Injection is pure gas
%[schedule.control.W.compi] = deal([1, 0]);
%% Set up initial state
ncomp = eos.fluid.getNumberOfComponents();
%7�����
z0 = [0.6, 0.3, 0.1, 0, 0, 0, 0];
T = 60 + 273.15;
p = 500*barsa;

pcswitch = 0;
state0 = initCompositionalState(G, p, T, [0.5, 0.5], z0, eos);
%% Simulate the schedule
pcswitch = 1;%PR-EOS����ë�����Ŀ���
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
        colorbar('vert') %����һ����ֱɫ��
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
%% %����states.pressure���ѷ�ѹ����
figure;
x=1:numel(states);
statespressure = [ ];
for step = 1:numel(states)
    statespressure(step,1) = states{step,1}.pressure(50,1);
end
plot(schstep,statespressure);
legend('fracturepressure');
%% %����wsgor��GOR��
figure;
x=1:numel(states);
wsgor = [ ] ;
for step = 1:numel(states)
    wsgor(step,1) = ws{step,1}.gor;
end
plot(schstep,wsgor);
legend('GOR');
%% %����states.rho�������ܶȡ�
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

%% --------------------------------------------------��ë����------------------------------------------------------
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm
global pcswitch fluidint;
%[filename,pathname]=uigetfile('.data','���data�ļ�');
%fn = fullfile(pathname,filename);%//�޸ģ�ѡ���ļ�
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

[schedulenon.control.W.components] = deal([1, 0, 0, 0, 0, 0, 0]);%7�����

%Injection is pure gas
%[schedule.control.W.compi] = deal([1, 0]);

%% Set up initial state
ncompnon = eosnon.fluid.getNumberOfComponents();
%7�����
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
%% %����states.x���ѷ��и���ֵ�Һ��ٷֱȡ�
figure;
compliquid = [ ] ;
for step = 1:numel(states)
    for i = 1:ncomp
    compliquid(step,i) = states{step,1}.x(50,i);
    end
end
plot(schstep,compliquid);
legend('Һ��', eos.fluid.names);
%% %����states.x�������и���ֵ�Һ��ٷֱȡ�
figure;
compliquidM = [ ] ;
for step = 1:numel(states)
    for i = 1:ncomp
    compliquidM(step,i) = states{step,1}.x(52,i);
    end
end
plot(schstep,compliquidM);
legend(eos.fluid.names);

%% %����statesnon.x���ѷ��и���ֵ�Һ��ٷֱȡ�
figure;
compliquidnon = [ ] ;
for step = 1:numel(statesnon)
    for i = 1:ncomp
    compliquidnon(step,i) = statesnon{step,1}.x(50,i);
    end
end
plot(schstep,compliquidnon);
legend('Һ��', eosnon.fluid.names);
%% %����statesnon.x�������и���ֵ�Һ��ٷֱȡ�
figure;
compliquidnonM = [ ] ;
for step = 1:numel(statesnon)
    for i = 1:ncomp
    compliquidnonM(step,i) = statesnon{step,1}.x(52,i);
    end
end
plot(schstep,compliquidnonM);
legend('Һ��', eosnon.fluid.names);
%% %����states.y���ѷ��и���ֵ�����ٷֱȡ�
figure;
compvapor = [ ] ;
for step = 1:numel(states)
    for i = 1:ncomp
    compvapor(step,i) = states{step,1}.y(50,i);
    end
end
plot(schstep,compvapor);
legend('����', eos.fluid.names);
%% %����states.y�������и���ֵ�����ٷֱȡ�
figure;
compvaporM = [ ] ;
for step = 1:numel(states)
    for i = 1:ncomp
    compvaporM(step,i) = states{step,1}.y(52,i);
    end
end
plot(schstep,compvaporM);
legend('����', eos.fluid.names);
%% %����statesnon.y���ѷ��и���ֵ�����ٷֱȡ�
figure;
compvapornon = [ ] ;
for step = 1:numel(statesnon)
    for i = 1:ncomp
    compvapornon(step,i) = statesnon{step,1}.y(50,i);
    end
end
plot(schstep,compvapornon);
legend('����', eosnon.fluid.names);
%% %����statesnon.y�������и���ֵ�����ٷֱȡ�
figure;
compvapornonM = [ ] ;
for step = 1:numel(statesnon)
    for i = 1:ncomp
    compvapornonM(step,i) = statesnon{step,1}.y(52,i);
    end
end
plot(schstep,compvapornonM);
legend('����', eosnon.fluid.names);
%% %����states.components���ѷ��и���ֵ�Ħ���ٷֱȡ�
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
%% %����states.components�������и���ֵ�Ħ���ٷֱȡ�
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
%% %����statesnon.components���ѷ��и���ֵ�Ħ���ٷֱȡ�
figure;
statescomponentsnon = [ ] ;
for step = 1:numel(statesnon)
    for i = 1:ncomp
    statescomponentsnon(step,i) = statesnon{step,1}.components(50,i);
    end
end
plot(schstep,statescomponentsnon);
legend(eosnon.fluid.names);
%% %����statesnon.components�������и���ֵ�Ħ���ٷֱȡ�
figure;
statescomponentsnonM = [ ] ;
for step = 1:numel(statesnon)
    for i = 1:ncomp
    statescomponentsnonM(step,i) = statesnon{step,1}.components(52,i);
    end
end
plot(schstep,statescomponentsnonM);
legend(eosnon.fluid.names);
%% %����states.pressure���ѷ�ѹ����
figure;
statespressurenon = [ ];
for step = 1:numel(statesnon)
    statespressurenon(step,1) = statesnon{step,1}.pressure(50,1);
end
plot(schstep,statespressurenon);
legend('fracturepressure');
%% %����wsgornon��GOR��
figure;
wsgornon = [ ] ;
for step = 1:numel(statesnon)
    wsgornon(step,1) = wsnon{step,1}.gor;
end
plot(schstep,wsgornon);
legend('GORnon');
%% %����statesnon.rho�������ܶȡ�
figure;
statesoilrhonon = [ ] ;
for step = 1:numel(statesnon)
    statesoilrhonon(step,1) = statesnon{step,1}.rho(50,1);
end
plot(schstep,statesoilrhonon);
legend('oilrhonon');
close all;

%% ����compareqOs_T
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
%% ����compareqGs_T
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
%% ����compareqgor
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
%% ����comparepressure
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
%% %����states.s�������еĺ��ͱ��Ͷȡ�
figure;
statesSO = [ ] ;
for step = 1:numel(states)
    statesSO(step,1) = states{step,1}.s(50,1);
end
plot(schstep,statesSO);
legend('SO');
%% %����statesnon.s���ѷ��еĺ��ͱ��Ͷȡ�
figure;
statesSOM = [ ] ;
for step = 1:numel(states)
    statesSOM(step,1) = states{step,1}.s(52,1);
end
plot(schstep,statesSOM);
legend('SO');
%% %����statesnon.s���ѷ��еĺ��ͱ��Ͷȡ�
figure;
statesSOnon = [ ] ;
for step = 1:numel(statesnon)
    statesSOnon(step,1) = statesnon{step,1}.s(50,1);
end
plot(schstep,statesSOnon);
legend('SO');
%% %����statesnon.s�������еĺ��ͱ��Ͷȡ�
figure;
statesSOnonM = [ ] ;
for step = 1:numel(statesnon)
    statesSOnonM(step,1) = statesnon{step,1}.s(52,1);
end
plot(schstep,statesSOnonM);
legend('SO');
%% %����states.Orho�������е������ܶȡ�
figure;
statesOrho = [ ] ;
for step = 1:numel(states)
    statesOrho(step,1) = states{step,1}.rho(50,1);
end
plot(schstep,statesOrho);
legend('Orho');
%% %����statesnon.Orho���ѷ��е������ܶȡ�
figure;
statesOrhoM = [ ] ;
for step = 1:numel(states)
    statesOrhoM(step,1) = states{step,1}.rho(52,1);
end
plot(schstep,statesOrhoM);
legend('Orho');
%% %����statesnon.Orho���ѷ��е������ܶȡ�
figure;
statesOrhonon = [ ] ;
for step = 1:numel(statesnon)
    statesOrhonon(step,1) = statesnon{step,1}.rho(50,1);
end
plot(schstep,statesOrhonon);
legend('Orho');
%% %����statesnon.Orho�������е������ܶȡ�
figure;
statesOrhononM = [ ] ;
for step = 1:numel(statesnon)
    statesOrhononM(step,1) = statesnon{step,1}.rho(52,1);
end
plot(schstep,statesOrhononM);
legend('Orho');
%% %����states.Grho�������е������ܶȡ�
figure;
statesGrho = [ ] ;
for step = 1:numel(states)
    statesGrho(step,1) = states{step,1}.rho(50,2);
end
plot(schstep,statesGrho);
legend('Grho');
%% %����statesnon.Grho���ѷ��е������ܶȡ�
figure;
statesGrhoM = [ ] ;
for step = 1:numel(states)
    statesGrhoM(step,1) = states{step,1}.rho(52,2);
end
plot(schstep,statesGrhoM);
legend('Grho');
%% %����statesnon.Grho���ѷ��е������ܶȡ�
figure; 
statesGrhonon = [ ] ;
for step = 1:numel(statesnon)
    statesGrhonon(step,1) = statesnon{step,1}.rho(50,2);
end
plot(schstep,statesGrhonon);
legend('Grho');
%% %����statesnon.Grho�������е������ܶȡ�
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