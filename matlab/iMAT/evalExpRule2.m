function [res] = evalExpRule2 (exp, x)
if (~( strcmp (exp , '') ) && ~isempty(exp))
    res = eval(exp);
else
    res = 0;
end
end

function [res] = and(a,b)
res = min(a,b);
end

function [res] = or(a,b)
res = max(a,b);

end



