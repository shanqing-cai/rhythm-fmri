function msglog(logFN, msg)
if ~isnumeric(logFN)
    logF = fopen(logFN, 'at');
else
    logF = logFN;
end

fprintf(logF, '%s\n', msg);

if ~isnumeric(logFN)
    fclose(logF);
end

fprintf(1, '%s\n', msg);

return