# Curves

**Curves** provides support for fitted curves (such as dose-response curves), 
including in-grid rendering, storing charts in cells, interactivity, and automatic fitting.

- Fitting: computing parameters of the specified function to best fit the data
    - Uses BFGS optimization algorithm (multi-threading for performance).
      For dose-response curves, we are typically fitting the sigmoid function
    - Ability to dynamically register custom fitting functions
    - Ability to get fitting performance characteristics (r-squared, classification, etc)
- Deep integration with the Datagrok grid
    - Either fitting on the fly, or using the supplied function + parameters
    - Multiple series in one cell
    - Ability to define chart, marker, or fitting options (such as fit function or marker color)
      on the column level, with the ability to override it on a grid cell or point level
    - Clicking a point in a chart within a grid makes it an outlier -> curve is re-fitted on the fly
    - Ability to specify a chart as "reference" so that it is shown on every other chart for comparison
- Ability to overlay curves from multiple grid cells (special viewer)
- Work with series stored in multiple formats (binary for performance, json for flexibility, etc)