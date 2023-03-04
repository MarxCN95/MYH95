function [schstep, qOs_D, qGs_D, qGs_T, qWs_D, qLs_D, fw, qOs_T, ...
    qWs_T, qLs_T, bhp] = GeneratePRO(ws, schedule, wellNow, varargin)
% 收集wellNow的生产动态数据
% SI单位制
opt = struct('Units', 'SI');
opt = merge_options(opt, varargin{:});
schstep    = schedule.step.val;
stepnum    = numel(ws);
qOs_D      = zeros(stepnum, 1);
qGs_D      = zeros(stepnum, 1);
qWs_D      = zeros(stepnum, 1);
bhp       = zeros(stepnum, 1);
for j = 1 : stepnum
    wL  = strcmp({ws{j, 1}.name}, wellNow);
    if  ~all(~wL)  % 有该口井
        if ws{j, 1}(wL).status % 该井开井
            qOs_D(j) = - ws{j, 1}(wL).qOs;  % 日产油
            bhp (j) =   ws{j, 1}(wL).bhp;  % 井底流压
            qGs_D(j)=  - ws{j, 1}(wL).qGs;  % 日产气
        end
    end
end
qLs_D         = qOs_D + qWs_D;             % 日产液
fw            = qWs_D./qLs_D;              % 含水率
fw(isnan(fw)) = 0;
qOs_T         = cumsum( qOs_D.* schstep ); % 累产油
qGs_T         = cumsum( qGs_D.* schstep ); % 累产气
qWs_T         = cumsum( qWs_D.* schstep ); % 累产水
qLs_T         = cumsum( qLs_D.* schstep ); % 累产液
% qOs_D = [0; qOs_D];
% qWs_D = [0; qWs_D];
% qLs_D = [0; qLs_D];
% fw    = [0;    fw];
% qOs_T = [0; qOs_T];
% qWs_T = [0; qWs_T];
% qLs_T = [0; qLs_T];
% qBHP  = [0;  qBHP];
schstep = cumsum(schstep);
if strcmp(opt.Units, 'Field')
    schstep = convertTo(schstep, day);
    qOs_D   = convertTo(qOs_D, meter^3/day);
    qGs_D   = convertTo(qGs_D, meter^3/day);
    qWs_D   = convertTo(qWs_D, meter^3/day);
    qLs_D   = convertTo(qLs_D, meter^3/day);
    bhp    = convertTo(bhp,  mega);
end

if strcmp(opt.Units, 'West')
    schstep = convertTo(schstep, day);
    qOs_D   = convertTo(qOs_D, stb/day);
    qGs_D   = convertTo(qGs_D, stb/day);
    qWs_D   = convertTo(qWs_D, stb/day);
    qLs_D   = convertTo(qLs_D, stb/day);
    bhp    = convertTo(bhp,  psia);
end