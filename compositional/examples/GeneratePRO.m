function [schstep, qOs_D, qGs_D, qGs_T, qWs_D, qLs_D, fw, qOs_T, ...
    qWs_T, qLs_T, bhp] = GeneratePRO(ws, schedule, wellNow, varargin)
% �ռ�wellNow��������̬����
% SI��λ��
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
    if  ~all(~wL)  % �иÿھ�
        if ws{j, 1}(wL).status % �þ�����
            qOs_D(j) = - ws{j, 1}(wL).qOs;  % �ղ���
            bhp (j) =   ws{j, 1}(wL).bhp;  % ������ѹ
            qGs_D(j)=  - ws{j, 1}(wL).qGs;  % �ղ���
        end
    end
end
qLs_D         = qOs_D + qWs_D;             % �ղ�Һ
fw            = qWs_D./qLs_D;              % ��ˮ��
fw(isnan(fw)) = 0;
qOs_T         = cumsum( qOs_D.* schstep ); % �۲���
qGs_T         = cumsum( qGs_D.* schstep ); % �۲���
qWs_T         = cumsum( qWs_D.* schstep ); % �۲�ˮ
qLs_T         = cumsum( qLs_D.* schstep ); % �۲�Һ
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