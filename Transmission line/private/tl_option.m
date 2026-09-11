function value = tl_option(rule,cfg,name,defaultValue,constraint)
%TL_OPTION Shared validation of scalar solver options; [] uses the default.
if isfield(cfg,name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
switch rule
    case 'text'
        if isstring(value) && isscalar(value), value = char(value); end
        if ~ischar(value) || size(value,1) ~= 1
            error('%s must be a text scalar.',name);
        end
        match = strcmpi(value,constraint);
        if ~any(match)
            error('%s must be one of: %s.',name,strjoin(constraint,', '));
        end
        value = constraint{find(match,1)};
        return;
    case 'logical'
        if ~islogical(value) || ~isscalar(value)
            error('%s must be a logical scalar.',name);
        end
        return;
end
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        isnan(value) || (isinf(value) && ~strcmp(rule,'positive-or-inf'))
    error('%s must be a real numeric scalar (finite unless Inf is allowed).',name);
end
switch rule
    case 'real'
        valid = true;
    case {'positive','positive-or-inf'}
        valid = value > 0;
    case 'nonnegative'
        valid = value >= 0;
    case 'integer'
        valid = value == round(value) && value >= constraint;
    case 'fraction'
        valid = value > 0 && value < 1;
    case 'fraction-or-zero'
        valid = value >= 0 && value < 1;
    otherwise
        error('Unknown option validation rule: %s.',rule);
end
if ~valid
    if strcmp(rule,'integer')
        error('%s must be an integer not smaller than %d.',name,constraint);
    end
    error('%s does not satisfy the %s constraint.',name,rule);
end
end
