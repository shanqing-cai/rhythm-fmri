function write_ost_to_file(ostMat, rmsSlopeWin, ostfn)
f = fopen(ostfn, 'wt');

fprintf(f, 'rmsSlopeWin = %.6f\n', rmsSlopeWin);
fprintf(f, '\n');

nRules = size(ostMat, 1);
fprintf(f, 'n = %d\n', nRules);
for i1 = 1 : nRules
    fprintf(f, '%d %d %.4f %.4f {}\n', ...
            ostMat(i1, 1), ostMat(i1, 2), ostMat(i1, 3), ostMat(i1, 4));
end

fclose(f);
return