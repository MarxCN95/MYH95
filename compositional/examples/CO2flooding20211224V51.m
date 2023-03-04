%% -------------------------------------------含有毛管力----------------------------------------------------------
%%  基础方案与只含毛管力的对比---2021年12月26日17:12:29
clear
mrstModule add ad-core ad-props mrst-gui compositional deckformat hfm linearsolvers
global pcswitch fluidint pcapillary;
global radius;
            %% ******************************------20----------*********

        radius = 20;
        
        includeWater = true;% Include aqueous phase
        %includeWater = false;%true;% Don't Include aqueous phase
        Shift = 0;  %critical properties shift
        
        %[filename,pathname]=uigetfile('.data','请打开data文件');
        %fn = fullfile(pathname,filename);%//修改，选择文件
        %fn    = fullfile('C:','Users','MYH','Desktop','YP2composition','YP2E300v1','YP2E300_E300.DATA');
        %fn    = fullfile('C:','Users','MYH','Desktop','YP2composition','YP2E300v2','YP2E300V2_E300.DATA');
        fn    = fullfile('C:','Users','MYH','Desktop','YP2composition','YP2E300v3fracturev3','YP2E300V3_E300.DATA');
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
        totTime = 5*year;
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
        %% Set up schedule
        dt = [0.001*86400;10*86400;30*86400;30*86400;30*86400;100*86400;100*86400;300*86400;300*86400;300*86400;300*86400;10*86400;30*86400;100*86400;100*86400];%100*86400;100*86400;100*86400;100*86400;100*86400;100*86400];%1*year;300*86400;300*86400;300*86400;300*86400;300*86400;300*86400;300*86400];%单位为秒，86400为一天。
        %      1           2         3        4       5        6         7 1year      8%       9         10        11         12    13       14        15
        %schedule = simpleSchedule(dt, 'W', W);
        dt = [];
        dt = [0.001*86400;1*86400;10*86400;30*86400;100*86400];
        for ti=6:9
            dt(ti) =1*year;
        end
        for ti=10:50
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
        
        %% start simulation
        pcswitch = 0;%PR-EOS考虑毛管力的开关
        pcapillary = 0 ;
        [wsb20, statesb20, repb20] = simulateScheduleAD(state0, model, schedule);
        pcswitch = 1;%PR-EOS考虑毛管力的开关
        pcapillary = 500 ;
        [wspcb20, statespcb20, reppcb20] = simulateScheduleAD(state0, model, schedule);

            %% ******************************------50----------*********

        radius = 50;
        
        includeWater = true;% Include aqueous phase
        %includeWater = false;%true;% Don't Include aqueous phase
        Shift = 0;  %critical properties shift
        
        %[filename,pathname]=uigetfile('.data','请打开data文件');
        %fn = fullfile(pathname,filename);%//修改，选择文件
        %fn    = fullfile('C:','Users','MYH','Desktop','YP2composition','YP2E300v1','YP2E300_E300.DATA');
        %fn    = fullfile('C:','Users','MYH','Desktop','YP2composition','YP2E300v2','YP2E300V2_E300.DATA');
        fn    = fullfile('C:','Users','MYH','Desktop','YP2composition','YP2E300v3fracturev3','YP2E300V3_E300.DATA');
        %fn    = fullfile('E:','Research','MRSTmodel','big grid3','1015_E300.DATA');
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
        totTime = 5*year;
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
        
        %% Set up schedule
        dt = [0.001*86400;10*86400;30*86400;30*86400;30*86400;100*86400;100*86400;300*86400;300*86400;300*86400;300*86400;10*86400;30*86400;100*86400;100*86400];%100*86400;100*86400;100*86400;100*86400;100*86400;100*86400];%1*year;300*86400;300*86400;300*86400;300*86400;300*86400;300*86400;300*86400];%单位为秒，86400为一天。
        %      1           2         3        4       5        6         7 1year      8%       9         10        11         12    13       14        15
        %schedule = simpleSchedule(dt, 'W', W);
        dt = [];
        dt = [0.001*86400;1*86400;10*86400;30*86400;100*86400];
        for ti=6:9
            dt(ti) =1*year;
        end
        for ti=10:50
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
        
        %% start simulation
        pcswitch = 0;%PR-EOS考虑毛管力的开关
        pcapillary = 0 ;
        [wsb50, statesb50, repb50] = simulateScheduleAD(state0, model, schedule);
        pcswitch = 1;%PR-EOS考虑毛管力的开关
        pcapillary = 500 ;
        [wspcb50, statespcb50, reppcb50] = simulateScheduleAD(state0, model, schedule);

            %% ******************************------0----------*********临界偏移为0， 仅毛管力100nm

        radius = 100;
        
        includeWater = true;% Include aqueous phase
        %includeWater = false;%true;% Don't Include aqueous phase
        Shift = 0;  %critical properties shift
        
        %[filename,pathname]=uigetfile('.data','请打开data文件');
        %fn = fullfile(pathname,filename);%//修改，选择文件
        %fn    = fullfile('C:','Users','MYH','Desktop','YP2composition','YP2E300v1','YP2E300_E300.DATA');
        %fn    = fullfile('C:','Users','MYH','Desktop','YP2composition','YP2E300v2','YP2E300V2_E300.DATA');
        fn    = fullfile('C:','Users','MYH','Desktop','YP2composition','YP2E300v3fracturev3','YP2E300V3_E300.DATA');
        %fn    = fullfile('E:','Research','MRSTmodel','big grid3','1015_E300.DATA');
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
        totTime = 5*year;
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
        
        %% Set up schedule
        dt = [0.001*86400;10*86400;30*86400;30*86400;30*86400;100*86400;100*86400;300*86400;300*86400;300*86400;300*86400;10*86400;30*86400;100*86400;100*86400];%100*86400;100*86400;100*86400;100*86400;100*86400;100*86400];%1*year;300*86400;300*86400;300*86400;300*86400;300*86400;300*86400;300*86400];%单位为秒，86400为一天。
        %      1           2         3        4       5        6         7 1year      8%       9         10        11         12    13       14        15
        %schedule = simpleSchedule(dt, 'W', W);
        dt = [];
        dt = [0.001*86400;1*86400;10*86400;30*86400;100*86400];
        for ti=6:9
            dt(ti) =1*year;
        end
        for ti=10:50
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
        
        %% start simulation
        pcswitch = 0;%PR-EOS考虑毛管力的开关
        pcapillary = 0 ;
        [wsb0, statesb0, repb0] = simulateScheduleAD(state0, model, schedule);
        pcswitch = 1;%PR-EOS考虑毛管力的开关
        pcapillary = 500 ;
        [wspcb100, statespcb100, reppcb100] = simulateScheduleAD(state0, model, schedule);
