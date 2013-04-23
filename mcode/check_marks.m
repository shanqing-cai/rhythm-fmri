function marksOkay = check_marks(this_utter, marks)
t_marks = nan(1, length(marks));

for i1 = 1 : numel(t_marks)
    t_marks(i1) = this_utter.(marks{i1});
end

marksOkay = isempty(find(diff(t_marks) <= 0));
return