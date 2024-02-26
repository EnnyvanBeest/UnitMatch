

function makepretty()
% set some graphical attributes of the current axis
set(gcf,'Color',[1 1 1])
set(get(gca, 'XLabel'), 'FontSize', 10);
set(get(gca, 'YLabel'), 'FontSize', 10);
set(gca, 'FontSize', 12,'TickDir','out');

set(get(gca, 'Title'), 'FontSize', 15);

ch = get(gca, 'Children');

for c = 1:length(ch)
    thisChild = ch(c);
    if strcmp('line', get(thisChild, 'Type')) || strcmp('errorbar', get(thisChild, 'Type'))
        if strcmp('.', get(thisChild, 'Marker'))
            set(thisChild, 'MarkerSize', 15);
        end
        if strcmp('-', get(thisChild, 'LineStyle'))
            set(thisChild, 'LineWidth', 2.0);
        end
        if strcmp('--', get(thisChild, 'LineStyle'))
            set(thisChild, 'LineWidth', 2.0);
        end       
    end
end
box off
