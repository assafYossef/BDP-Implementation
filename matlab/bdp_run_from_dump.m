function [e, info, data] = bdp_run_from_dump(dumpDir, opts)
%BDP_RUN_FROM_DUMP Run the solver on the matrices written by `bdp --dump <dir>`.
%
%   Lets you develop / debug the MATLAB side without the C++ process. The C++
%   program writes every matrix it would push through the MATLAB Engine as a
%   plain text file <dumpDir>/<name>.txt. The result is written back to
%   <dumpDir>/e_opt.txt so that `bdp <mesh> --load-metric <dumpDir>/e_opt.txt`
%   can visualise it.

if nargin < 2, opts = bdp_default_options(); end
names = {'F','F2E','EV','e0','cA','cV','cW','wF','isBnd','Vc','K0','Eref','isCone','coneAngle'};
for i = 1:numel(names)
    fn = fullfile(dumpDir, [names{i} '.txt']);
    if exist(fn, 'file')
        info = dir(fn);
        if info.bytes == 0
            data.(names{i}) = zeros(0, 1);      % e.g. no constrained vertices
        else
            data.(names{i}) = load(fn, '-ascii');
        end
    elseif ~any(strcmp(names{i}, {'Eref','isCone','coneAngle'}))
        data.(names{i}) = zeros(0, 1);  % e.g. no constrained vertices
    end
end
optFile = fullfile(dumpDir, 'variant.txt');
if exist(optFile, 'file') && ~isfield(opts, 'variantFromCaller')
    fid = fopen(optFile, 'r'); v = strtrim(fgetl(fid)); fclose(fid);
    opts.variant = v;
end
[e, info] = bdp_solve(data, opts);
fid = fopen(fullfile(dumpDir, 'e_opt.txt'), 'w');
fprintf(fid, '%.17g\n', e);
fclose(fid);
end
