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


## Software requirement

The `.fig` files should be opened with MATLAB. Compatibility with GNU Octave is not guaranteed.
