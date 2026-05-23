function plot_branch_json(jsonfile, yField, linewidthStable, linewidthUnstable, colorStable, colorUnstable)
%PLOT_BRANCH_JSON Plot a BifurcationKit branch saved as JSON.
%
% Example:
%   figure; hold on;
%   plot_branch_json('data/br_primary.json', 'meanu');
%   plot_branch_json('data/br_from_first_bp.json', 'meanu', 2, 1, [0 0.45 0.74], [0.85 0.33 0.1]);
%   hold off;

if nargin < 2 || isempty(yField)
    yField = 'meanu';
end
if nargin < 3 || isempty(linewidthStable)
    linewidthStable = 3;
end
if nargin < 4 || isempty(linewidthUnstable)
    linewidthUnstable = 1;
end
if nargin < 5 || isempty(colorStable)
    colorStable = [0 0 0];
end
if nargin < 6 || isempty(colorUnstable)
    colorUnstable = [1 0 0];
end

S = jsondecode(fileread(jsonfile));
br = S.branch;
n = numel(br);

x = zeros(n,1);
y = zeros(n,1);
stable = false(n,1);

for i = 1:n
    x(i) = br(i).param;
    y(i) = br(i).(yField);

    if isfield(br(i), 'stable')
        stable(i) = br(i).stable ~= 0;
    else
        stable(i) = false;
    end
end

hold on;

for i = 1:n-1
    if stable(i) && stable(i+1)
        plot(x(i:i+1), y(i:i+1), '-', ...
            'Color', colorStable, ...
            'LineWidth', linewidthStable);
    else
        plot(x(i:i+1), y(i:i+1), '-', ...
            'Color', colorUnstable, ...
            'LineWidth', linewidthUnstable);
    end
end

xlabel('\omega_{d0}');
ylabel(yField, 'Interpreter', 'none');
grid on;
box on;

end
