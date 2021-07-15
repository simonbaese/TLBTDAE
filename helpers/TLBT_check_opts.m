function tlbt_opts = TLBT_check_opts(tlbt_opts)
%% Checks opts fields for TLBT.
%
% Inputs:
% tlbt_opts    option struct used for TLBT model order reduction
%
% Author: Simon Bäse (2021)
%
% See LICENSE.md

if not(isfield(tlbt_opts,'tl'))
    tlbt_opts.tl = 0;
else
    if not(isnumeric(tlbt_opts.tl)) && not(islogical(tlbt_opts.tl))
        error('tlbt_opts.tl parameter must be logical or numeric.');
    end
end

if not(isfield(tlbt_opts,'tf'))
    warning('tlbt_opts.tf not set. Fallback to traditional balanced truncation');
    tlbt_opts.tl = 0;
else
    if not(isnumeric(tlbt_opts.tf)) || not(tlbt_opts.tf >= 0)
        error('tlbt_opts.tf parameter must be numeric and positive.');
    end
end

if not(isfield(tlbt_opts,'adiplot'))
    tlbt_opts.adiplot = 0;
else
    if not(isnumeric(tlbt_opts.adiplot)) && not(islogical(tlbt_opts.adiplot))
        error('tlbt_opts.adiplot parameter must be logical or numeric.');
    end
end

if not(isfield(tlbt_opts,'analyze'))
    tlbt_opts.analyze = 0;
else
    if not(isnumeric(tlbt_opts.analyze)) && not(islogical(tlbt_opts.analyze))
        error('tlbt_opts.analyze parameter must be logical or numeric.');
    end
end

if not(isfield(tlbt_opts,'simulate'))
    tlbt_opts.simulate = 0;
else
    if not(isnumeric(tlbt_opts.simulate)) && not(islogical(tlbt_opts.simulate))
        error('tlbt_opts.simulate parameter must be logical or numeric.');
    end
end

if tlbt_opts.simulate
    if not(isfield(tlbt_opts,'tf'))
        error('tlbt_opts.tf required for simulation.');
    end     
    if not(isfield(tlbt_opts,'horizon'))
        if not(isfield(tlbt_opts,'tf'))
            tlbt_opts.horizon = 5;
        else
            tlbt_opts.horizon = tlbt_opts.tf;
        end
    else
        if not(isnumeric(tlbt_opts.horizon)) || not(tlbt_opts.horizon >= tlbt_opts.tf)
            error('tlbt_opts.horizon parameter must be numeric and larger than tlbt_opts.tf.');
        end
    end
    if not(isfield(tlbt_opts,'input'))
        tlbt_opts.input = 'impulse';
    else
        if not(strcmp(tlbt_opts.input, 'impulse')) && not(strcmp(tlbt_opts.input, 'step'))
            warning('tlbt_opts.input must be impulse or step. Falling back to impulse.');
            tlbt_opts.input = 'impulse';
        end
    end
end

end

