function write_ost_to_file(ostMat, rmsSlopeWin, maxInterval_5_10, ostfn)
f = fopen(ostfn, 'wt');

fprintf(f, 'rmsSlopeWin = %.6f\n', rmsSlopeWin);
fprintf(f, '\n');

nRules = size(ostMat, 1);
fprintf(f, 'n = %d\n', nRules);
for i1 = 1 : nRules
    fprintf(f, '%d %d %.4f %.4f {}\n', ...
            ostMat(i1, 1), ostMat(i1, 2), ostMat(i1, 3), ostMat(i1, 4));
end

fprintf(f, '\n');
fprintf(f, 'n = 1\n');
fprintf(f, '5 %.6f 10\n', maxInterval_5_10);

fclose(f);
return