# Review: 11-compute (Scientific computing)

Fresh-eyes review (2026-11-22, ~22:45-23:20) of the 7 delivered clips: feature, descriptions.md block, and contact
sheets at the key lines and results of each delivered take (every delivered video matched its latest take before
the review). Judged on: is the one idea clear to someone new to Datagrok, is each line true when it shows (data
claims checked against the frame), plain wording, jumps.

| NN | Clip | Verdict |
|---|---|---|
| 01 | Run a model in Diff Studio (0:45) | OK. The hub shows the templates and the library; Prey 30 shows in the table's first row; alpha's tooltip reads "The maximum prey per capita growth rate"; Finish 150 gives 1,501 rows. |
| 02 | Write your own model (1:05) | OK. After Refresh the logistic curve rises from 10 toward K = 1000; with r at 1.5 it levels off by t = 8. |
| 03 | Fit a model to data (1:22) | refilmed (the coordinator's question). The Fit view is a table view, so with Browse open the toolbox took the left dock: the form was squeezed and the result grid's charts cut at about 27 of 36 h. The status bar's Toolbox toggle only handed the dock back to Browse, so the clip now closes Browse from its title bar right after the model opens ("Close Browse for more room"), and the Fit view comes up without a side panel: the charts show all three doses with the model over the data, and the next output's chart starts to show. The best fit comes first (loss 2.39E-2, then 3.08E-2, 6.75E-2), clearance 15.82. The "Target: what the model should match" line shows a moment before the Target block scrolls into view (the next click brings it). |
| 04 | Find which inputs matter with sensitivity analysis (1:12) | refilmed, the same fix. Its claims, from the new run (Sobol is random): clearance has the longest bars (first order 0.74, total 0.78, on Central conc.), the correlation tooltip reads Pearson R -0.760 for clearance against the central level. The form's labels are still cut ("Absorption ra..."), less than before. The feature's maintainer note said 100 samples; it is 1000, fixed. |
| 05 | Run a model from Model Hub (0:40) | title card fixed (re-rendered from the same take): "Model Hub lists the models published on your server" → "the models you can run on your server" (the cards are Compute2's catalog of functions tagged as models, here the ones the installed plugins bring); the description says the same. Max distance 12.61 and max height 3.962 after the changes. |
| 06 | Save a model to your files (0:40) | OK. Predator-prey.ivp is in My files and opens the model again (the view's tab reads "Lotka-Volterra (2)": the model's own name). The opened view brings the toolbox back on the left. |
| 07 | Export a model's equations to LaTeX (0:46) | OK. With Standalone off the source starts at \textbf{Lotka-Volterra} instead of \documentclass; the file name is Lotka-Volterra-2026-11-22.tex. |

Verdicts: 4 OK, 2 refilmed (03, 04), 1 title card fixed (05).

## Data claims

All true: 01 (Prey 30 in the first row, 1,501 rows after Finish 150), 02 (the curve levels off at K), 03 (best fit
first, the model drawn over three doses of data), 04 (clearance first in both Sobol charts, R -0.76), 05 (max
distance and height change with the inputs), 06 (the file is in My files), 07 (the preamble goes with Standalone
off).

## Titles and descriptions

Plain and specific; each names where the thing is (Browse > Apps > Diff Studio or Model Hub). Changed: 05's title
card and description ("published" was loose), 03/04 m:ss.

## Systematic

- Diff Studio's Fit and Sensitivity views are table views: opened while Browse is docked, they bring the toolbox
  into the same dock. To give a view the full width on camera, close Browse with the x on its title bar (hover the
  panel first: its title bar icons show only under the pointer while filming). The sidebar's compass, clicked with
  the pointer over the sidebar, does not close a docked Browse (`simple_mode.dart`: the stripe's `hover` class
  sends it to `hideSlidingToolbox`), and the status bar's Toolbox toggle only swaps the dock back to Browse.
- 06 and 07 open the model from Browse as well; 06's reopened model shows the toolbox, harmless there.
