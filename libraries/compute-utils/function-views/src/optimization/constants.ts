// Specific optimization constants

/** Dock ratios */
export enum DOCK_RATIO {
    PC_PLOT = 0.5,
    GRAPH = 0.6,
}

/** Size constants */
export const ROW_HEIGHT = 25;

/** */
enum HELP_ITEMS {
  GOAL = '`Goal`',
  SAMPLES = '`Samples`',
};

/** Starting help markdown */
export const STARTING_HELP = `# Optimization

Function [optimization](https://en.wikipedia.org/wiki/Mathematical_optimization) involves
determining the input values that either maximize or minimize a given output.

1. Use switchers to specify **target** scalar outputs:   
   * If you select one of them, then the [Nelder-Mead](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method)
method runs. Set type of optimization in ${HELP_ITEMS.GOAL}   
   * If you select several scalars, then the [Monte Carlo](https://en.wikipedia.org/wiki/Monte_Carlo_method)
method performs

2. Apply switchers to set inputs for optimization. Once you choose any of them, set a range of its variation.

3. Press **Run** or <i class="fas fa-play"></i> on the top panel.

4. Press <i class="grok-icon fal fa-question"></i> on the top panel to learn more about 
[Optimization](https://datagrok.ai/help/compute.md#optimization).

## Single target optimization

If you select one **target** output, then the 
[Nelder-Mead](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method) method method runs. 
Press <i class="grok-icon fal fa-cog"></i> to edit its settings.

You will get

**TO ADD**

## Multi target analysis

If you select more than one **target** output, then the 
[Monte Carlo](https://en.wikipedia.org/wiki/Monte_Carlo_method) method runs. 
It computates the function multiple times with varying inputs, 
and analyzes the relationship between them and outputs. Press <i class="grok-icon fal fa-cog"></i>
to set a number of the function evaluations.

You will get

  * [PC plot](https://datagrok.ai/help/visualize/viewers/pc-plot) visualizing multivariate data
and providing variations of the selected inputs & outputs

  * [Line chart](https://datagrok.ai/help/visualize/viewers/line-chart) or 
[Scatterplot](https://datagrok.ai/help/visualize/viewers/scatter-plot) (dependently on the varied
inputs count) showing a behavior of each output separately

  * [Grid](https://datagrok.ai/help/visualize/viewers/grid) containing inputs and outputs values 
of each function evaluation`;
