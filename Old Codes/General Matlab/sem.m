% Aaron Hossack
% Standard Error of the Mean
function[y] = sem(x)
if size(x, 1) == 1
    x = x'; % make into a column
end
y = std(x) ./ sqrt(size(x, 1));