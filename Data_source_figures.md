# MATLAB figure source files (`.fig`)

This link  https://drive.google.com/drive/folders/1D4_CwcYtc263bpFLvDXxWxADq71PmtLW?usp=sharing  contains the MATLAB `.fig` files corresponding to the figures used in the manuscript and in the supplementary material.

The files are provided to make the figures reproducible and to allow readers to inspect and, when possible, extract the numerical data stored inside the MATLAB figure objects.

## Naming convention

The figure-source files follow the numbering used in the manuscript:

```text
Fig1.fig
Fig2.fig
Fig3a.fig
Fig3b.fig
...
FigS1.fig
FigS1a.fig
...
```

For multipanel figures, the panel letter is appended to the figure number, for example:

```text
Fig6a.fig
Fig6b.fig
Fig6c.fig
```

Supplementary figures are labelled with the prefix `FigS`.

## Opening a figure in MATLAB

A `.fig` file can be opened directly in MATLAB with:

```matlab
openfig('Fig1.fig');
```

To open it without displaying the window:

```matlab
fig = openfig('Fig1.fig','invisible');
```

After opening the file, the graphical objects can be inspected using `findall` or `findobj`.

## Extracting line data

For a simple figure containing curves plotted with `plot`, the data are usually stored in objects of type `Line`.

```matlab
fig = openfig('Fig1.fig','invisible');

ax = findall(fig,'Type','axes');
lines = findall(ax,'Type','line');

for i = 1:numel(lines)
    x = lines(i).XData;
    y = lines(i).YData;

    T = table(x(:), y(:), 'VariableNames', {'x','y'});
    writetable(T, sprintf('Fig1_line_%02d.csv', i));
end

close(fig);
```

This creates one CSV file for each curve stored in the MATLAB figure.

## Extracting scatter data

For figures created with `scatter`, the data are usually stored in objects of type `Scatter`.

```matlab
fig = openfig('Fig1.fig','invisible');

ax = findall(fig,'Type','axes');
sc = findall(ax,'Type','scatter');

for i = 1:numel(sc)
    x = sc(i).XData;
    y = sc(i).YData;

    T = table(x(:), y(:), 'VariableNames', {'x','y'});
    writetable(T, sprintf('Fig1_scatter_%02d.csv', i));
end

close(fig);
```

## Extracting surface or image data

For figures created with commands such as `surf`, `mesh`, `pcolor`, `imagesc`, or similar routines, the data may be stored as `Surface` or `Image` objects.

### Surface objects

```matlab
fig = openfig('Fig1.fig','invisible');

ax = findall(fig,'Type','axes');
surfaces = findall(ax,'Type','surface');

for i = 1:numel(surfaces)
    X = surfaces(i).XData;
    Y = surfaces(i).YData;
    Z = surfaces(i).ZData;
    C = surfaces(i).CData;

    save(sprintf('Fig1_surface_%02d.mat', i), 'X', 'Y', 'Z', 'C');
end

close(fig);
```

### Image objects

```matlab
fig = openfig('Fig1.fig','invisible');

ax = findall(fig,'Type','axes');
images = findall(ax,'Type','image');

for i = 1:numel(images)
    X = images(i).XData;
    Y = images(i).YData;
    C = images(i).CData;

    save(sprintf('Fig1_image_%02d.mat', i), 'X', 'Y', 'C');
end

close(fig);
```

For image-like plots, `CData` contains the matrix displayed in the figure. If the figure was saved from an already rasterized image, `CData` may contain only pixel values, not the original numerical simulation data.

## General MATLAB extraction script

The following MATLAB function extracts the most common data types from a `.fig` file and saves them in a dedicated output folder.

Save it as `extract_fig_data.m`.

```matlab
function extract_fig_data(figfile, outdir)
%EXTRACT_FIG_DATA Extract common data objects from a MATLAB .fig file.
%
% Usage:
%   extract_fig_data('Fig1.fig', 'extracted_data')
%
% The function extracts:
%   - Line objects       -> CSV files with x, y, and possibly z columns
%   - Scatter objects    -> CSV files with x, y, and possibly z columns
%   - Surface objects    -> MAT files containing X, Y, Z, C
%   - Image objects      -> MAT files containing X, Y, C
%
% Notes:
%   The availability of numerical data depends on how the original figure was created.
%   If the figure contains only a raster image, the original raw data may not be recoverable.

    if nargin < 2 || isempty(outdir)
        outdir = 'extracted_data';
    end

    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    [~, basename, ~] = fileparts(figfile);

    fig = openfig(figfile, 'invisible');
    cleanupObj = onCleanup(@() close(fig));

    axes_list = findall(fig, 'Type', 'axes');

    % Remove legend axes if present
    keep = true(size(axes_list));
    for k = 1:numel(axes_list)
        tag = get(axes_list(k), 'Tag');
        if strcmpi(tag, 'legend')
            keep(k) = false;
        end
    end
    axes_list = axes_list(keep);

    for a = 1:numel(axes_list)
        ax = axes_list(a);

        % Lines
        lines = findall(ax, 'Type', 'line');
        for i = 1:numel(lines)
            x = lines(i).XData(:);
            y = lines(i).YData(:);

            try
                z = lines(i).ZData(:);
            catch
                z = [];
            end

            if isempty(z)
                T = table(x, y, 'VariableNames', {'x','y'});
            else
                T = table(x, y, z, 'VariableNames', {'x','y','z'});
            end

            fname = fullfile(outdir, sprintf('%s_axis%02d_line%02d.csv', basename, a, i));
            writetable(T, fname);
        end

        % Scatter
        scatters = findall(ax, 'Type', 'scatter');
        for i = 1:numel(scatters)
            x = scatters(i).XData(:);
            y = scatters(i).YData(:);

            try
                z = scatters(i).ZData(:);
            catch
                z = [];
            end

            if isempty(z)
                T = table(x, y, 'VariableNames', {'x','y'});
            else
                T = table(x, y, z, 'VariableNames', {'x','y','z'});
            end

            fname = fullfile(outdir, sprintf('%s_axis%02d_scatter%02d.csv', basename, a, i));
            writetable(T, fname);
        end

        % Surfaces
        surfaces = findall(ax, 'Type', 'surface');
        for i = 1:numel(surfaces)
            X = surfaces(i).XData;
            Y = surfaces(i).YData;
            Z = surfaces(i).ZData;
            C = surfaces(i).CData;

            fname = fullfile(outdir, sprintf('%s_axis%02d_surface%02d.mat', basename, a, i));
            save(fname, 'X', 'Y', 'Z', 'C');
        end

        % Images
        images = findall(ax, 'Type', 'image');
        for i = 1:numel(images)
            X = images(i).XData;
            Y = images(i).YData;
            C = images(i).CData;

            fname = fullfile(outdir, sprintf('%s_axis%02d_image%02d.mat', basename, a, i));
            save(fname, 'X', 'Y', 'C');
        end
    end

    fprintf('Data extracted from %s into folder: %s\n', figfile, outdir);
end
```

Example:

```matlab
extract_fig_data('Fig1.fig', 'extracted_data');
extract_fig_data('FigS1.fig', 'extracted_data');
```

To extract all `.fig` files in the current folder:

```matlab
files = dir('*.fig');

for k = 1:numel(files)
    extract_fig_data(files(k).name, 'extracted_data');
end
```

## Important limitations

The `.fig` files preserve the MATLAB graphics objects used to create the figures. Therefore, the possibility of recovering numerical data depends on the content of the figure:

- if the figure contains curves produced with `plot`, the vectors are usually stored as `XData` and `YData`;
- if the figure contains scatter plots, the points are usually stored as `XData` and `YData`;
- if the figure contains surfaces or heat maps, the arrays may be stored as `XData`, `YData`, `ZData`, and `CData`;
- if the figure contains a raster image imported into MATLAB, only the displayed pixel matrix may be available, not the original simulation data.

For complete reproducibility, the `.fig` files should be used together with the numerical scripts and, when available, the raw simulation data.

## GitHub notes

MATLAB `.fig` files are binary files. GitHub can store them, but it cannot show meaningful text diffs.

If the files are large, it is recommended to track them with Git LFS:

```bash
git lfs install
git lfs track "*.fig"
git add .gitattributes
git add *.fig
git commit -m "Add MATLAB figure source files"
```

A recommended repository structure is:

```text
figures/
  Fig1.fig
  Fig2.fig
  Fig3a.fig
  Fig3b.fig
  FigS1.fig

scripts/
  extract_fig_data.m

README_figures.md
```

## Software requirement

The `.fig` files should be opened with MATLAB. Compatibility with GNU Octave is not guaranteed.
