function [res] = evalExpRule (exp, x)
if (~( strcmp (exp , '') ) && ~isempty(exp))
    res = eval(exp);
else
    res = 0;
end
end

function [res] = and(a,b)
%     res = min(a,b);
%     if (a==0 & b==1) | (a==1 & b==0)
%         res = 1;
%     end

if (a==1 & b==0) | (a==0 & b==1) | (a==-1 & b==0) | (a==0 & b==-1) | (a==1 & b==-1) | (a==-1 & b==1)
    res = 0;
end
if (a==1 & b==1)
    res = 1;
end
if (a==-1 & b==-1)
    res=-1;
end
if (a==0 & b==0)
    res=0;
end
end

function [res] = or(a,b)%didn't take care of 1 | -1
res = max(a,b);
if (a==0 & b==-1) | (a==-1 & b==0)
    res = -1;
end
if (a==1 & b==-1)| (a==-1 & b==1)
    res = 0;
end
end



