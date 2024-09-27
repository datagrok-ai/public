// Specific sensitivity analysis constants

/** Dock ratios */
export enum DOCK_RATIO {
    COR_PLOT = 0.5,
    PC_PLOT = 0.6,
    GRAPH = PC_PLOT,
    BAR_CHART = 0.2,
  }

/** Size constants */
export const ROW_HEIGHT = 25;

enum HELP_ITEMS {
  METHOD = '`Method`',
  SAMPLES = '`Samples`',
};

/** Starting help markdown */
export const STARTING_HELP = `# Sensitivity Analysis

Sensitivity Analysis runs the computation multiple times with varying inputs, 
and analyzes the relationship between inputs and outputs.

Use one of the following methods:

* [Monte Carlo](https://datagrok.ai/help/compute/function-analysis#monte-carlo) explores a function 
at randomly taken points

* [Sobol](https://datagrok.ai/help/compute/function-analysis#sobol)
decomposes output variance into fractions, which can be attributed to inputs

* [Grid](https://datagrok.ai/help/compute/function-analysis#grid) studies a function at the points of a grid with the specified step

## Monte Carlo

Once you've chosen it in ${HELP_ITEMS.METHOD}

* Set in ${HELP_ITEMS.SAMPLES} the number of random points 

* Use switchers to specify varied inputs and outputs to be analyzed

* Press **Run** or <i class="fas fa-play"></i> on the top panel. You will get 

  * [Correlation plot](https://datagrok.ai/help/visualize/viewers/correlation-plot) for exploring
correlations between varied inputs and the specified outputs

  * [PC plot](https://datagrok.ai/help/visualize/viewers/pc-plot) visualizing multivariate data
and providing variations of the selected inputs & outputs

  * [Line chart](https://datagrok.ai/help/visualize/viewers/line-chart) or 
[Scatterplot](https://datagrok.ai/help/visualize/viewers/scatter-plot) (dependently on the varied
inputs count) showing a behavior of each output separately

  * [Grid](https://datagrok.ai/help/visualize/viewers/grid) containing inputs and outputs values 
of each function evaluation

* Press <i class="grok-icon fal fa-question"></i> on the top panel to learn more about 
[Sensitivity Analysis](https://datagrok.ai/help/compute/function-analysis#sensitivity-analysis)

## Sobol

This method performs 
[variance-based sensitivity analysis](https://en.wikipedia.org/wiki/Variance-based_sensitivity_analysis). 
It provides the same visualizations 
as **Monte Carlo** and [bar charts](https://datagrok.ai/help/visualize/viewers/bar-chart) showing 
Sobol' indices:

* First-order indices indicate the contribution to the output variance of varying each input alone

* Total-order indices measure the contribution to the output variance of each input,
 including all variance caused by its interactions with any other inputs
 
 ## Grid
 
 This method provides the same visualization as the previous ones. Unlike them, it evaluates a function 
 at the points of uniform grid within the specified ranges.`;
