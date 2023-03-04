function [stable, x, y] = phaseStabilityTest(eos, z, p, T, x, y, state, cells,varargin)
% Perform a phase stability test for a mixture
%
% SYNOPSIS:
%   [stable, x, y] = phaseStabilityTest(eos, z, p, T, x, y)
%
% DESCRIPTION:
%   Perform a phase stability test for a multicomponent 
%
% PARAMETERS:
%   eos - EquationOfState derived class instance.
%   z   - Composition as a matrix with number of rows equal to the number
%         of components.
%   p   - Pressures as a column vector
%   T   - Temperatures as a column vector
%   x   - Initial guess for liquid mole fractions as a matrix with number
%         of rows equal to the number of components. For the most general
%         case, this is equal to z.
%   y   - Initial guess for vapor mole fractions as a matrix with number
%         of rows equal to the number of components. For the most general
%         case, this is equal to z.
%
% OPTIONAL PARAMETERS:
%
%  'tol_equil'   - Tolerance for equilibrium constants in stability test.
%                  Uses the natural logarithm of K_i squared as a norm.
%
%  'tol_trivial' - Tolerance for the trivial test
%
%  'solve_both'  - Perform a stability test for a second liquid-like phase
%                  for cells already determined to be unstable by the
%                  second vapor-like phase. This has slightly higher
%                  computational cost since the stability must be checked
%                  twice, but the phase fraction estimates are more
%                  accurate. Enabled by default if x, y are requested as
%                  outputs.
%
% RETURNS:
%   stable - Indicator. If the entry for cell i is true, then the mixture
%            is stable (single-phase) for the p, T regime in that cell.
%   x      - Estimates of liquid mole fractions. Equal to z for stable
%            cells.
%   y      - Estimates of vapor mole fractions. Equal to z for stable
%            cells.
%
% SEE ALSO:
%   `EquationOfStateModel`

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
    opt = struct('tol_equil',       1e-10, ...
                 'alpha',           [], ...
                 'tol_trivial',     1e-5, ...
                 'solve_both',      nargout > 1);
    assert(all(p > 0 & isfinite(p)), 'Positive, finite pressures required for phase stability test.');
    if numel(varargin)
        opt = merge_options(opt, varargin{:});
    end
    if nargin < 5
        x = z;
        y = z;
    end
    S_tol = opt.tol_trivial;
    nc = numel(p);
    active = true(nc, 1);

    z = ensureMinimumFraction(z);
    global pcswitch
    if pcswitch == 1
    [A_ij, Bi] = eos.getMixingParameters(p, T, eos.fluid.acentricFactors, false);
    [Si_L, A_L, B_L] = eos.getPhaseMixCoefficients(x, A_ij, Bi);
    [Si_V, A_V, B_V] = eos.getPhaseMixCoefficients(y, A_ij, Bi);
    Z_L = eos.computeCompressibilityZ(p, x, A_L, B_L, Si_L, Bi, 1);
    Z_V = eos.computeCompressibilityZ(p, y, A_V, B_V, Si_V, Bi, 0);
    rho_L = eos.PropertyModel.computeDensity(p, x, Z_L, T, 1);
rho_V = eos.PropertyModel.computeDensity(p, y, Z_V, T, 0);
[Pcap] = computePcapfinal(eos.fluid, rho_L,rho_V,x,y);
Pcap = Pcap*0;
%Pcap = Pcap(cells, :);
p_V = p + Pcap;%�޸ģ� IFT����ë����
    [y, S_V, isTrivialV] = checkStability_V(eos, z, y, p_V, T, true, active, opt, state);%�޸�
    else
         [y, S_V, isTrivialV] = checkStability(eos, z, y, p, T, true, active, opt);%�޸�   
    end
    V_stable = S_V <= 1.1 + S_tol | isTrivialV;%���Ŵ�0.01%
    
    if ~opt.solve_both
        active = V_stable;
    end
    [x, S_L, isTrivialL] = checkStability(eos, z, x, p, T, false, active, opt);
    L_stable = S_L <= 1 + S_tol | isTrivialL;

    stable = (L_stable & V_stable);
    x(stable, :) = z(stable, :);
    y(stable, :) = z(stable, :);
    
    bad = ~isfinite(S_L);
    if any(bad)
        for i = 1:numel(x)
            x(bad, :) = z(bad, :);
        end
    end
    bad = ~isfinite(S_V);
    if any(bad)
        for i = 1:numel(y)
            y(bad, :) = z(bad, :);
        end
    end
end

function [xy, S, trivialSolution, K] = checkStability(eos, z, xy, p, T, insidePhaseIsVapor, active, opt)
    tol_trivial = opt.tol_trivial;
    tol_equil = opt.tol_equil;
    nc = numel(p);

    [A_ij, Bi] = eos.getMixingParameters(p, T, eos.fluid.acentricFactors, false);
    
    f_z = getFugacity(eos, A_ij, Bi, z, p, insidePhaseIsVapor);
    
    K = estimateEquilibriumWilson(eos, p, T);

        
    % xy is either x or y, depending on context for phase test
    S = zeros(nc, 1);
    trivialSolution = false(nc, 1);
    if ~any(active)
        return
    end
    for stepNo = 1:20000
        zi = z(active, :);
        ki = K(active, :);
        A_ij_loc = cellfun(@(x) x(active, :), A_ij, 'UniformOutput', false);
        Bi_loc = Bi(active, :);
        
        if insidePhaseIsVapor
            xy_loc = ki.*zi;
        else
            xy_loc = zi./ki;
        end
        S_loc = sum(xy_loc, 2);
        xy_loc = bsxfun(@rdivide, xy_loc, S_loc);
        
%         if xy_loc(1,7)==1.892192800794860e-04
%           warning('numelnan have NaN'); 
%     end                   % 2021��11��4��22:25:06  �޸�
%     xy_loc(1,7)==1.892192800794860e-04�����
        f_xy = getFugacity(eos, A_ij_loc, Bi_loc, xy_loc, p(active), ~insidePhaseIsVapor);
        f_zi = f_z(active, :);
        if insidePhaseIsVapor
            % f_xy is vapor fugacity
            R = bsxfun(@rdivide, f_zi./f_xy, S_loc);
        else
            % f_xy is liquid fugacity
            R = bsxfun(@times, f_xy./f_zi, S_loc);
        end
    numelnan = numel(find(isnan(R)));
    if numelnan > 1
       warning('numelnan have NaN'); 
    end     % 2021��11��4��22:25:06  �޸� ������� xy_loc��������޴��ֵ������R����NaN��ֵ
        K(active, :) = K(active, :).*R;
        R_norm = sum((R-1).^2, 2);
        K_norm = sum(log(K(active, :)).^2, 2);

        trivial = K_norm < tol_trivial;
        converged = R_norm < tol_equil;
        done = trivial | converged;
        S(active) = S_loc;
        trivialSolution(active) = trivial;
        xy(active, :) = xy_loc;
        active(active) = ~done;
        if all(done)
            dispif(mrstVerbose() > 1, 'Stability done in %d iterations\n', stepNo);
            break
        end
    end
    if ~all(done)
        warning('Stability test did not converge');
    end
end


function f = getFugacity(model, A_ij, Bi, xy, p, isLiquid)
    [Si, A, B] = model.getPhaseMixCoefficients(xy, A_ij, Bi);
    Z = model.computeCompressibilityZ(p, xy, A, B, Si, Bi, isLiquid);
    f = model.computeFugacity(p, xy, Z, A, B, Si, Bi);
end
%% ����Ϊ�޸�
function [xy, S, trivialSolution, K] = checkStability_L(eos, z, xy, p, T, insidePhaseIsVapor, active, opt)
    tol_trivial = opt.tol_trivial;
    tol_equil = opt.tol_equil;
    nc = numel(p);

    [A_ij, Bi] = eos.getMixingParameters(p, T, eos.fluid.acentricFactors, false);
    
    f_z = getFugacity(eos, A_ij, Bi, z, p, insidePhaseIsVapor);
    
    K = estimateEquilibriumWilson(eos, p, T);
    % xy is either x or y, depending on context for phase test
    S = zeros(nc, 1);
    trivialSolution = false(nc, 1);
    if ~any(active)
        return
    end
    for stepNo = 1:20000
        zi = z(active, :);
        ki = K(active, :);
        A_ij_loc = cellfun(@(x) x(active, :), A_ij, 'UniformOutput', false);
        Bi_loc = Bi(active, :);
        
        if insidePhaseIsVapor
            xy_loc = ki.*zi;
        else
            xy_loc = zi./ki;
        end
        S_loc = sum(xy_loc, 2);
        xy_loc = bsxfun(@rdivide, xy_loc, S_loc);
        
        f_xy = getFugacity(eos, A_ij_loc, Bi_loc, xy_loc, p(active), ~insidePhaseIsVapor);
        f_zi = f_z(active, :);
        if insidePhaseIsVapor
            % f_xy is vapor fugacity
            R = bsxfun(@rdivide, f_zi./f_xy, S_loc);
        else
            % f_xy is liquid fugacity
            R = bsxfun(@times, f_xy./f_zi, S_loc);
        end

        K(active, :) = K(active, :).*R;
        R_norm = sum((R-1).^2, 2);
        K_norm = sum(log(K(active, :)).^2, 2);

        trivial = K_norm < tol_trivial;
        converged = R_norm < tol_equil;
        done = trivial | converged;
        S(active) = S_loc;
        trivialSolution(active) = trivial;
        xy(active, :) = xy_loc;
        active(active) = ~done;
        if all(done)
            dispif(mrstVerbose() > 1, 'Stability done in %d iterations\n', stepNo);
            break
        end
    end
    if ~all(done)
        warning('Stability test did not converge');
    end
end

%% ����Ϊ�޸�2020-11-12
function [xy, S, trivialSolution, K] = checkStability_V(eos, z, xy, p, T, insidePhaseIsVapor, active, opt, state)
    tol_trivial = opt.tol_trivial;
    tol_equil = opt.tol_equil;
    nc = numel(p);
%     global pcswitch  pcapillary;
%     if pcswitch == 1
%         p_V = p + pcapillary;
%     end %//���ë����
    [A_ij, Bi] = eos.getMixingParameters(p, T, eos.fluid.acentricFactors, false);
    
    f_z = getFugacity(eos, A_ij, Bi, z, p, insidePhaseIsVapor);
    
    K = estimateEquilibriumWilson(eos, p, T);
    % xy is either x or y, depending on context for phase test
    S = zeros(nc, 1);
    trivialSolution = false(nc, 1);
    if ~any(active)
        return
    end
    for stepNo = 1:20000
        zi = z(active, :);
        ki = K(active, :);
        A_ij_loc = cellfun(@(x) x(active, :), A_ij, 'UniformOutput', false);
        Bi_loc = Bi(active, :);
        
        if insidePhaseIsVapor
            xy_loc = ki.*zi;
        else
            xy_loc = zi./ki;
        end
        S_loc = sum(xy_loc, 2);
        xy_loc = bsxfun(@rdivide, xy_loc, S_loc);
        
        f_xy = getFugacity_V(eos, A_ij_loc, Bi_loc, xy_loc, p(active), ~insidePhaseIsVapor, state);
        f_zi = f_z(active, :);
        if insidePhaseIsVapor
            % f_xy is vapor fugacity
            R = bsxfun(@rdivide, f_zi./f_xy, S_loc);
        else
            % f_xy is liquid fugacity
            R = bsxfun(@times, f_xy./f_zi, S_loc);
        end

        K(active, :) = K(active, :).*R;
        R_norm = sum((R-1).^2, 2);
        K_norm = sum(log(K(active, :)).^2, 2);

        trivial = K_norm < tol_trivial;
        converged = R_norm < tol_equil;
        done = trivial | converged;
        S(active) = S_loc;
        trivialSolution(active) = trivial;
        xy(active, :) = xy_loc;
        active(active) = ~done;
        if all(done)
            dispif(mrstVerbose() > 1, 'Stability done in %d iterations\n', stepNo);
            break
        end
    end
    if ~all(done)
        warning('Stability test did not converge');
    end
end


function f = getFugacity_V(model, A_ij, Bi, y, p, isLiquid, state)            
    [Si, A, B] = model.getPhaseMixCoefficients(y, A_ij, Bi);
    Z = model.computeCompressibilityZ(p, y, A, B, Si, Bi, isLiquid);
    f = model.computeFugacity(p, y, Z, A, B, Si, Bi);
end