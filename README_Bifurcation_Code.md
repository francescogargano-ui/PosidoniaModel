# Posidonia 1D bifurcation code

Clean Julia code for the 1D bifurcation analysis of the reduced *Posidonia oceanica* model.

## Files

```text
posidonia_bifurcation_clean.jl
plot_branch_json.m
```

## Run the Julia code

```bash
julia posidonia_bifurcation_clean.jl
```

The script computes the primary homogeneous branch, detects the first bifurcation point, continues one branch emerging from it, and saves the results as JSON files:

```text
data/br_primary.json
data/br_from_first_bp.json
```

## Plot branches in MATLAB

The branches saved in JSON can be plotted with:

```matlab
figure; hold on;

plot_branch_json('data/br_primary.json', 'meanu');
plot_branch_json('data/br_from_first_bp.json', 'meanu', ...
                 2, 1, [0 0.447 0.741], [0.85 0.325 0.098]);

xlabel('\omega_{d0}');
ylabel('mean(n)');
title('Bifurcation diagram');

hold off;
```

Other fields saved in the JSON can also be used, for example:

```matlab
plot_branch_json('data/br_primary.json', 'maxn');
plot_branch_json('data/br_primary.json', 'amp');
```
