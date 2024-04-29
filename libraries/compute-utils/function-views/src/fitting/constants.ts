// Specific optimization constants

/** Dock ratios */
export enum GRID_SIZE {
  LOSS_COL_WIDTH = 60,
  ROW_HEIGHT = 180,
  LOSS_GRAPH_WIDTH = 360,
  GOF_VIEWER_WIDTH = 240,
}

/** UI titles */
export enum TITLE {
   ITER = 'Iteration',
   LOSS = 'Loss',
   LOSS_GRAPH = 'Fitting profile',
   TARGET = 'Target',
   FIT = 'Fit',
   METHOD = 'method',
   VALUE = 'Value',
   OBJECTIVE = TARGET,
   SAMPLES = 'samples',
   ID = 'id',
 };

/** Items for the fitting help */
enum HELP_ITEMS {
  GOAL = '`Goal`',
  FIT = '`' + TITLE.FIT + '`',
  MIN = '`min`',
  MAX = '`max`',
  TARGET = '`' + TITLE.TARGET + '`',
  METHOD = '`' + TITLE.METHOD + '`',
  CONTEXT = '`Context Panel (F4)`',
  SAMPLES = '`' + TITLE.SAMPLES + '`',
};

/** Starting help markdown */
export const STARTING_HELP = `## Fitting

Use fitting to solve an inverse problem: find input conditions leading to specified output constraints. 
It computes inputs minimizing losses measured by 
[root mean square deviation](https://en.wikipedia.org/wiki/Root-mean-square_deviation).

1. In the ${HELP_ITEMS.FIT} block, use switchers to specify inputs to be found:
   * Set ${HELP_ITEMS.MIN} and ${HELP_ITEMS.MAX} values for each selected item. They define the variation range
   * Set values of all other inputs

2. Set output constraints in the ${HELP_ITEMS.TARGET} block:
   * Use switchers to specify target outputs
   * Set target value for each selected item

3. Choose ${HELP_ITEMS.METHOD}. Press <i class="grok-icon fal fa-cog"></i>
to specify its settings.

4. Specify the number of inputs sets to be obtained (in the ${HELP_ITEMS.SAMPLES} field).

5. Press **Run** or <i class="fas fa-play"></i> on the top panel to perform fitting. You will get:   

   * [Grid](https://datagrok.ai/help/visualize/viewers/grid) containing obtained inputs and loss values

   * [PC plot](https://datagrok.ai/help/visualize/viewers/pc-plot) providing variations of the fitted inputs
and obtained loss

   * [Line chart](https://datagrok.ai/help/visualize/viewers/line-chart) showing loss minimization

   * viewers illustrating goodness of fit

6. Click on the grid row, and explore the selected set:

   * The line chart shows the loss function behavior. It illustrates the process of the current sample fitting

   * The goodness of fit viewer visualizes deviation of the obtained output from its target value

   * Open ${HELP_ITEMS.CONTEXT}. You will get the function run corresponding to the selected row

8. Press <i class="grok-icon fal fa-question"></i> on the top panel to learn more about 
[parameters optimization](https://datagrok.ai/help/compute/#input-parameter-optimization).

## Learn more

* [Parameters optimization](https://datagrok.ai/help/compute/#input-parameter-optimization)
* [Nelder-Mead method](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method)`;
