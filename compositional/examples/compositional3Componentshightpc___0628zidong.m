%% -------------------------------------------含有毛管力----------------------------------------------------------
clear
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm
global pcswitch fluidint;
%[filename,pathname]=uigetfile('.data','请打开data文件');
%fn = fullfile(pathname,filename);%//修改，选择文件
fn    = fullfile('C:','Users','MYH','Desktop','big grid3','1015_E300.DATA');
%fn    = fullfile('C:','Users','myh','Desktop','0905','0623_E300.DATA');
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
[schedule.control.W.components] = deal([1, 0, 0, 0, 0, 0, 0]);%7个组分

%Injection is pure gas
%[schedule.control.W.compi] = deal([1, 0]);
%% Set up initial state
ncomp = eos.fluid.getNumberOfComponents();
%7个组分
z0 = [0.333, 0.333, 0.333, 0, 0, 0, 0];
T = 60 + 273.15;
p = 60*barsa;

pcswitch = 0;
state0 = initCompositionalState(G, p, T, [0.5, 0.5], z0, eos);
%% Simulate the schedule
pcswitch = 1;%PR-EOS考虑毛管力的开关
[ws, states, rep] = simulateScheduleAD(state0, model, schedule);
%%
[schstep, qOs_D, qGs_D, qGs_T, qWs_D, qLs_D, fw, qOs_T, ...
    qWs_T, qLs_T, bhp] = GeneratePRO(ws, schedule, 'PROD', 'Units', 'Field');

%%
%% --------------------------------------------------无毛管力------------------------------------------------------
%% --------------------------------------------------无毛管力------------------------------------------------------
%% --------------------------------------------------无毛管力------------------------------------------------------
%% --------------------------------------------------无毛管力------------------------------------------------------
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm
global pcswitch fluidint;
%[filename,pathname]=uigetfile('.data','请打开data文件');
%fn = fullfile(pathname,filename);%//修改，选择文件
%fnnon    = fullfile('C:','Users','myh','Desktop','0905','0623_E300.DATA');
fnnon = fn;
decknon = readEclipseDeck(fnnon);
decknon = convertDeckUnits(decknon);

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
pcswitch = 0;
state0non = initCompositionalState(Gnon, p, T, [0.5, 0.5], z0, eosnon);
%% Simulate the schedule
pcswitch = 0;%%%%%%%%%%%%%%%考虑毛管力的开关%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[wsnon, statesnon, repnon] = simulateScheduleAD(state0non, modelnon, schedulenon);
%%
[schstep, qOs_Dnon, qGs_Dnon, qGs_Tnon, qWs_Dnon, qLs_Dnon, fwnon, qOs_Tnon, ...
    qWs_Tnon, qLs_Tnon, bhpnon] = GeneratePRO(wsnon, schedulenon, 'PROD', 'Units', 'Field');
endtime = datestr(now);
save('60bar');
%% -------------------------------------------含有毛管力----------------------------------------------------------
clear
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm
global pcswitch fluidint;
%[filename,pathname]=uigetfile('.data','请打开data文件');
%fn = fullfile(pathname,filename);%//修改，选择文件
fn    = fullfile('C:','Users','myh','Desktop','big grid3','1015_E300.DATA');
%fn    = fullfile('C:','Users','myh','Desktop','0905','0623_E300.DATA');
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
[schedule.control.W.components] = deal([1, 0, 0, 0, 0, 0, 0]);%7个组分

%Injection is pure gas
%[schedule.control.W.compi] = deal([1, 0]);
%% Set up initial state
ncomp = eos.fluid.getNumberOfComponents();
%7个组分
z0 = [0.333, 0.333, 0.333, 0, 0, 0, 0];
T = 60 + 273.15;
p = 100*barsa;

pcswitch = 0;
state0 = initCompositionalState(G, p, T, [0.5, 0.5], z0, eos);
%% Simulate the schedule
pcswitch = 1;%PR-EOS考虑毛管力的开关
[ws, states, rep] = simulateScheduleAD(state0, model, schedule);
%%
[schstep, qOs_D, qGs_D, qGs_T, qWs_D, qLs_D, fw, qOs_T, ...
    qWs_T, qLs_T, bhp] = GeneratePRO(ws, schedule, 'PROD', 'Units', 'Field');

%%
%% --------------------------------------------------无毛管力------------------------------------------------------
%% --------------------------------------------------无毛管力------------------------------------------------------
%% --------------------------------------------------无毛管力------------------------------------------------------
%% --------------------------------------------------无毛管力------------------------------------------------------
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm
global pcswitch fluidint;
%[filename,pathname]=uigetfile('.data','请打开data文件');
%fn = fullfile(pathname,filename);%//修改，选择文件
%fnnon    = fullfile('C:','Users','myh','Desktop','0905','0623_E300.DATA');
fnnon = fn;
decknon = readEclipseDeck(fnnon);
decknon = convertDeckUnits(decknon);

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
pcswitch = 0;
state0non = initCompositionalState(Gnon, p, T, [0.5, 0.5], z0, eosnon);
%% Simulate the schedule
pcswitch = 0;%%%%%%%%%%%%%%%考虑毛管力的开关%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[wsnon, statesnon, repnon] = simulateScheduleAD(state0non, modelnon, schedulenon);
%%
[schstep, qOs_Dnon, qGs_Dnon, qGs_Tnon, qWs_Dnon, qLs_Dnon, fwnon, qOs_Tnon, ...
    qWs_Tnon, qLs_Tnon, bhpnon] = GeneratePRO(wsnon, schedulenon, 'PROD', 'Units', 'Field');
endtime = datestr(now);
save('100bar');