function print_progress(fraction_done)
    persistent t
    if fraction_done >= 1 % Task completed
        t = [];
        fprintf('\n');
        return
    end
    if isempty(t)
        fprintf('%5.2f%%', 0);
        t = 1;
    end
    fprintf('\b\b\b\b\b\b%5.2f%%', 100*fraction_done);
end
