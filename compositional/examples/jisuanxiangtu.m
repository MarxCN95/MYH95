%% -------------------------------------------����ë����----------------------------------------------------------
clear
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm
global pcswitch fluidint;
%[filename,pathname]=uigetfile('.data','���data�ļ�');
%fn = fullfile(pathname,filename);%//�޸ģ�ѡ���ļ�
fn    = fullfile('C:','Users','myh','Desktop','0829','0623_E300.DATA');

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
z0 = [0.333, 0.333, 0.333, 0, 0, 0, 0];
T = 60 + 273.15;
p = 70*barsa;

pcswitch = 0;
state0 = initCompositionalState(G, p, T, [0.5, 0.5], z0, eos);
%% Simulate the schedule
pcswitch = 1;%PR-EOS����ë�����Ŀ���
eosmodel = model.EOSModel;
[stable, x, y] = phaseStabilityTest(eosmodel, z0, p, T, z0, z0, state0);
            