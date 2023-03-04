%% -------------------------------------------����ë����----------------------------------------------------------
clear
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm
global pcswitch fluidint;
%[filename,pathname]=uigetfile('.data','���data�ļ�');
%fn = fullfile(pathname,filename);%//�޸ģ�ѡ���ļ�
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
[schedule.control.W.components] = deal([1, 0, 0, 0, 0, 0, 0]);%7�����

%Injection is pure gas
%[schedule.control.W.compi] = deal([1, 0]);
%% Set up initial state
ncomp = eos.fluid.getNumberOfComponents();
%7�����
z0 = [0.333, 0.333, 0.333, 0, 0, 0, 0];
T = 60 + 273.15;
p = 60*barsa;

pcswitch = 0;
state0 = initCompositionalState(G, p, T, [0.5, 0.5], z0, eos);
%% Simulate the schedule
pcswitch = 1;%PR-EOS����ë�����Ŀ���
[ws, states, rep] = simulateScheduleAD(state0, model, schedule);
%%
[schstep, qOs_D, qGs_D, qGs_T, qWs_D, qLs_D, fw, qOs_T, ...
    qWs_T, qLs_T, bhp] = GeneratePRO(ws, schedule, 'PROD', 'Units', 'Field');

%%
%% --------------------------------------------------��ë����------------------------------------------------------
%% --------------------------------------------------��ë����------------------------------------------------------
%% --------------------------------------------------��ë����------------------------------------------------------
%% --------------------------------------------------��ë����------------------------------------------------------
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm
global pcswitch fluidint;
%[filename,pathname]=uigetfile('.data','���data�ļ�');
%fn = fullfile(pathname,filename);%//�޸ģ�ѡ���ļ�
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

[schedulenon.control.W.components] = deal([1, 0, 0, 0, 0, 0, 0]);%7�����

%Injection is pure gas
%[schedule.control.W.compi] = deal([1, 0]);

%% Set up initial state
ncompnon = eosnon.fluid.getNumberOfComponents();
%7�����
pcswitch = 0;
state0non = initCompositionalState(Gnon, p, T, [0.5, 0.5], z0, eosnon);
%% Simulate the schedule
pcswitch = 0;%%%%%%%%%%%%%%%����ë�����Ŀ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[wsnon, statesnon, repnon] = simulateScheduleAD(state0non, modelnon, schedulenon);
%%
[schstep, qOs_Dnon, qGs_Dnon, qGs_Tnon, qWs_Dnon, qLs_Dnon, fwnon, qOs_Tnon, ...
    qWs_Tnon, qLs_Tnon, bhpnon] = GeneratePRO(wsnon, schedulenon, 'PROD', 'Units', 'Field');
endtime = datestr(now);
save('60bar');
%% -------------------------------------------����ë����----------------------------------------------------------
clear
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm
global pcswitch fluidint;
%[filename,pathname]=uigetfile('.data','���data�ļ�');
%fn = fullfile(pathname,filename);%//�޸ģ�ѡ���ļ�
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
[schedule.control.W.components] = deal([1, 0, 0, 0, 0, 0, 0]);%7�����

%Injection is pure gas
%[schedule.control.W.compi] = deal([1, 0]);
%% Set up initial state
ncomp = eos.fluid.getNumberOfComponents();
%7�����
z0 = [0.333, 0.333, 0.333, 0, 0, 0, 0];
T = 60 + 273.15;
p = 100*barsa;

pcswitch = 0;
state0 = initCompositionalState(G, p, T, [0.5, 0.5], z0, eos);
%% Simulate the schedule
pcswitch = 1;%PR-EOS����ë�����Ŀ���
[ws, states, rep] = simulateScheduleAD(state0, model, schedule);
%%
[schstep, qOs_D, qGs_D, qGs_T, qWs_D, qLs_D, fw, qOs_T, ...
    qWs_T, qLs_T, bhp] = GeneratePRO(ws, schedule, 'PROD', 'Units', 'Field');

%%
%% --------------------------------------------------��ë����------------------------------------------------------
%% --------------------------------------------------��ë����------------------------------------------------------
%% --------------------------------------------------��ë����------------------------------------------------------
%% --------------------------------------------------��ë����------------------------------------------------------
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm
global pcswitch fluidint;
%[filename,pathname]=uigetfile('.data','���data�ļ�');
%fn = fullfile(pathname,filename);%//�޸ģ�ѡ���ļ�
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

[schedulenon.control.W.components] = deal([1, 0, 0, 0, 0, 0, 0]);%7�����

%Injection is pure gas
%[schedule.control.W.compi] = deal([1, 0]);

%% Set up initial state
ncompnon = eosnon.fluid.getNumberOfComponents();
%7�����
pcswitch = 0;
state0non = initCompositionalState(Gnon, p, T, [0.5, 0.5], z0, eosnon);
%% Simulate the schedule
pcswitch = 0;%%%%%%%%%%%%%%%����ë�����Ŀ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[wsnon, statesnon, repnon] = simulateScheduleAD(state0non, modelnon, schedulenon);
%%
[schstep, qOs_Dnon, qGs_Dnon, qGs_Tnon, qWs_Dnon, qLs_Dnon, fwnon, qOs_Tnon, ...
    qWs_Tnon, qLs_Tnon, bhpnon] = GeneratePRO(wsnon, schedulenon, 'PROD', 'Units', 'Field');
endtime = datestr(now);
save('100bar');