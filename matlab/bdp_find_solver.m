function [name, reason] = bdp_find_solver(request, verbose)
%BDP_FIND_SOLVER Pick the conic solver for the SCP subproblems.
%
%   [name, reason] = bdp_find_solver(request)
%
%   request : 'auto' | 'mosek' | 'cvx' | 'coneprog' | 'newton'
%   name    : the solver that will be used
%   reason  : one line explaining the choice (printed by bdp_solve)
%
%   'auto' tries, in this order:
%     1. MOSEK called directly (mosekopt). mosekopt is often NOT on the MATLAB path
%        even when MOSEK is installed: CVX ships its own copy of MOSEK in
%        <cvx>/mosek/<arch> and only puts it on the path while cvx_end is running.
%        We therefore also look there, and in the standard MOSEK installation
%        folders. A tiny LP is solved to make sure the licence allows direct use
%        (CVX's bundled MOSEK may be licensed for use through CVX only).
%     2. CVX, which uses whatever solver you selected with cvx_solver / cvx_setup
%        (e.g. MOSEK) - exactly how the reference script solves its subproblems.
%     3. coneprog (MATLAB Optimization Toolbox R2020b+).
%     4. the built-in Newton barrier solver.

if nargin < 2, verbose = false; end
request = lower(request);
isOctave = exist('OCTAVE_VERSION', 'builtin') > 0;
reason = '';

switch request
    case 'newton'
        name = 'newton'; reason = 'requested'; return;
    case 'mosek'
        [ok, why] = mosek_usable();
        if ~ok, error('bdp:solver', 'MOSEK cannot be used directly: %s', why); end
        name = 'mosek'; reason = why; return;
    case 'cvx'
        if ~cvx_available(), error('bdp:solver', 'CVX not found on the MATLAB path (run cvx_setup)'); end
        name = 'cvx'; reason = sprintf('CVX with its default solver (%s)', cvx_default_solver()); return;
    case 'coneprog'
        if isOctave || exist('coneprog', 'file') == 0
            error('bdp:solver', 'coneprog (Optimization Toolbox R2020b+) not available');
        end
        name = 'coneprog'; reason = 'requested'; return;
    case 'auto'
        [ok, why] = mosek_usable();
        if ok, name = 'mosek'; reason = why; return; end
        msg = why;
        if cvx_available()
            name = 'cvx';
            reason = sprintf('direct MOSEK not usable (%s); using CVX with its default solver (%s)', ...
                             msg, cvx_default_solver());
            return;
        end
        if ~isOctave && exist('coneprog', 'file') > 0
            name = 'coneprog';
            reason = sprintf('direct MOSEK not usable (%s) and CVX not found', msg);
            return;
        end
        name = 'newton';
        reason = sprintf('no MOSEK, CVX or coneprog found (%s)', msg);
    otherwise
        error('bdp:solver', 'Unknown solver "%s" (auto | mosek | cvx | coneprog | newton)', request);
end
if verbose, fprintf('[bdp] solver: %s - %s\n', name, reason); end
end

% =============================================================================
function [ok, why] = mosek_usable()
ok = false;
where = '';
if exist('mosekopt', 'file') == 0
    dirs = candidate_mosek_dirs();
    for i = 1:numel(dirs)
        if ~isempty(dir(fullfile(dirs{i}, 'mosekopt.*')))
            addpath(dirs{i});                       % keep it for this MATLAB session
            where = dirs{i};
            break;
        end
    end
    if exist('mosekopt', 'file') == 0
        why = 'mosekopt not found on the path, in the CVX folder or in the MOSEK install folders';
        return;
    end
else
    where = fileparts(which('mosekopt'));
end
% licence / sanity check: minimise x s.t. x >= 1
try
    prob.c = 1; prob.a = sparse(1, 1, 1); prob.blc = 1; prob.buc = inf; prob.blx = 0; prob.bux = inf;
    [rc, res] = mosekopt('minimize echo(0)', prob);
    ok = (rc == 0) && abs(res.sol.itr.xx - 1) < 1e-6;
    if ok
        why = sprintf('mosekopt found in %s', where);
    else
        why = sprintf('mosekopt in %s returned code %d (licence?)', where, rc);
    end
catch err
    why = sprintf('mosekopt in %s failed: %s', where, err.message);
end
if ~ok && ~isempty(where) && ~isempty(strfind(path, where))
    % do not leave an unusable copy on the path (CVX manages its own copy)
    if exist('cvx_where', 'file') > 0 && ~isempty(strfind(where, cvx_where()))
        rmpath(where);
    end
end
end

function dirs = candidate_mosek_dirs()
dirs = {};
% 1) MOSEK bundled with CVX:  <cvx>/mosek/<arch>
if exist('cvx_where', 'file') > 0
    root = cvx_where();
    sub = dir(fullfile(root, 'mosek'));
    for i = 1:numel(sub)
        if sub(i).isdir && sub(i).name(1) ~= '.'
            dirs{end+1} = fullfile(root, 'mosek', sub(i).name); %#ok<AGROW>
        end
    end
    dirs{end+1} = fullfile(root, 'mosek');
end
% 2) standalone MOSEK installations: <base>/<version>/toolbox/<release>
bases = {};
if ispc
    bases = {fullfile(getenv('ProgramFiles'), 'Mosek'), 'C:\Program Files\Mosek', fullfile(getenv('USERPROFILE'), 'mosek')};
else
    bases = {fullfile(getenv('HOME'), 'mosek'), '/opt/mosek', '/usr/local/mosek'};
end
for b = 1:numel(bases)
    vers = dir(bases{b});
    for v = numel(vers):-1:1                         % newest version first
        if ~vers(v).isdir || vers(v).name(1) == '.', continue; end
        tb = fullfile(bases{b}, vers(v).name, 'toolbox');
        rel = dir(tb);
        for r = numel(rel):-1:1
            if rel(r).isdir && rel(r).name(1) ~= '.'
                dirs{end+1} = fullfile(tb, rel(r).name); %#ok<AGROW>
            end
        end
    end
end
end

function ok = cvx_available()
ok = exist('cvx_begin', 'file') > 0 && exist('cvx_end', 'file') > 0;
end

function s = cvx_default_solver()
s = 'unknown';
try
    s = cvx_solver();                                % current CVX solver name
catch
end
end
