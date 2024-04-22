// Specific optimization constants

/** Dock ratios */
export enum DOCK_RATIO {
    PC_PLOT = 0.5,
    LOSS_PLOT = 0.5,
}

/** Size constants */
export const ROW_HEIGHT = 25;

/** Items for the fitting help */
enum HELP_ITEMS {
  GOAL = '`Goal`',
  SAMPLES = '`Samples`',
  FIT = '`Fit`',
  MIN = '`min`',
  MAX = '`max`',
  TO_GET = '`Target`',
  METHOD = '`method`',
  CONTEXT = '`Context Panel (F4)`',
};

/** Starting help markdown */
export const STARTING_HELP = `## Fitting

Use fitting to solve an inverse problem: find input conditions leading to specified output constraints. 
It computes inputs minimizing deviation measured by [loss function](https://en.wikipedia.org/wiki/Loss_function).

1. In the ${HELP_ITEMS.FIT} block, use switchers to specify inputs to be found:
   * Set ${HELP_ITEMS.MIN} and ${HELP_ITEMS.MAX} values for each selected item. They define the variation range
   * Set values of all other inputs

2. Set output constraints in the ${HELP_ITEMS.TO_GET} block:
   * Use switchers to specify target outputs
   * Set target value for each selected item

3. Choose ${HELP_ITEMS.METHOD}. Press <i class="grok-icon fal fa-cog"></i>
to edit its settings.

4. Press **Run** or <i class="fas fa-play"></i> on the top panel to perform fitting. You will get:  

   * [Line chart](https://datagrok.ai/help/visualize/viewers/line-chart) showing loss function minimization

   * [PC plot](https://datagrok.ai/help/visualize/viewers/pc-plot) providing variations of the fitted inputs
and loss function

   * [Grid](https://datagrok.ai/help/visualize/viewers/grid) containing values of the fitted inputs and
the loss function for each iteration

5. Explore viewers vizualizing the fitting process

6. To get an evaluation of particular interest: 
   * Click on grid row
   * Open ${HELP_ITEMS.CONTEXT}. You will get the function run corresponding to the selected row

7. Press <i class="grok-icon fal fa-question"></i> on the top panel to learn more about 
[parameters optimization](https://datagrok.ai/help/compute/#input-parameter-optimization).

## Learn more

* [Parameters optimization](https://datagrok.ai/help/compute/#input-parameter-optimization)
* [Nelder-Mead method](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method)
* [Gradient descent](https://en.wikipedia.org/wiki/Gradient_descent)`;

export enum TITLE {
  ITER = 'Iteration',
  LOSS = 'Loss',
};
