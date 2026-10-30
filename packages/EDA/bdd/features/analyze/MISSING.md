# EDA analyze features — what is still missing

One section per group of a round. Each lists the scenarios that need a step the library lacks (as Gherkin,
with the signal it would read), the questions for the help or the product, and the findings.

## EDA dialogs

Features: `reduce-dimensionality.feature`, `clustering.feature`, `hypothesis-testing.feature`
(30 Oct 2026, local stand).

### Needs a step (EDA dialogs)

- **T-test: why Run is disabled.** The dialog binds the reason to Run as a tooltip (`ui.tooltip.bind` in
  `ttest/ttest-ui.ts` `updateRunButtonState`: "Column has 3 groups (3+). Use ANOVA.", "Each group needs at
  least 2 observations…", "Alpha must be strictly between 0 and 1", the Student/unequal-variances advice),
  and the red icon next to Category carries it the same way. `{element} should be disabled because {string}`
  reads `aria-description` or `title`, which neither has. Either the dialog sets the reason as the button's
  `aria-description` (then the existing step reads it), or a step hovers the button and reads the tooltip:

  ```gherkin
  Scenario: A category with three groups says why Run is disabled
    When user picks "ML > Analyze > Group Comparison > T-test..." from the top menu
    And user selects "site" in Category input in "Two-sample t-test" dialog
    Then Run button in "Two-sample t-test" dialog should be disabled because "Column has 3 groups (3+). Use ANOVA."
    And no errors should have been logged
  ```

- **Use WebGPU** (Reduce Dimensionality's UMAP settings, Markov clustering). Headless Chromium has no WebGPU
  adapter, so the input is off and disabled ("WebGPU is not available"). A run with an adapter needs a
  capability gate like the service gates, `the browser has a WebGPU adapter`
  (`navigator.gpu?.requestAdapter()`), before the scenario that switches it on.

### Questions (EDA dialogs)

- **ML > Cluster > DBSCAN... and MCL... have no help page.** `explore/cluster-data.md` describes only k-means
  (Cluster Data...); DBSCAN's parameters are explained only as a postprocessing of Reduce Dimensionality
  (`explore/dim-reduction.md`). The features claim current behaviour. Should both items get a section?
- **DBSCAN dialog captions are raw parameter names**: "Df", "X Col", "Y Col" (`dbScan` in EDA
  `src/package.ts` declares no caption for them), next to the captioned "Epsilon" and "Minimum points".
  Intended, or should they read Table / X / Y?
- **Reduce Dimensionality, Postprocessing**: the stand offers "DBSCAN clustering" and "None"; the help names
  them "DBSCAN" and "None" and adds "Radial Coloring", which the NLP package registers (not installed on the
  local stand). Should the help say Radial Coloring needs NLP? The help also calls the dialog "Dimensionality
  Reduction parameter dialog"; its title is "Reduce Dimensionality".
- **Markov clustering: OK with no column picked** is enabled (unlike Reduce Dimensionality, whose editor it
  extends). What should it do? Not claimed.

### Findings 30 Oct 2026 (EDA dialogs)

- **Reduce Dimensionality: OK with no column is greyed out but not disabled.** Steps: open cars; ML > Reduce
  Dimensionality...; leave Columns empty. Measured: the OK button has the CSS classes `enabled disabled`
  (greyed, `pointer-events: none`), but `button.disabled` is false and it has no `aria-disabled`, so assistive
  technology and the BDD state reading see an enabled button (the run reads "disabled expected … false"; a
  probe confirmed `disabled=false aria=null`). Expected: OK disabled, as the T-test dialog does with
  `runBtn.disabled = true`. Source: EDA `src/package.ts:220`, `okButton.classList.add('disabled')` in
  `reduceDimensionality`'s `validate`. The help says nothing either way; "OK disabled with no input" is the
  platform's dialog convention. `reduce-dimensionality.feature`, scenario "With no column picked OK is
  disabled", `@known-failure` on that step (runs p1d, p1e).
- **Markov clustering on number columns only never finishes.** Steps: open iris; ML > Cluster > MCL...;
  Columns: tick Sepal.Length, Sepal.Width, Petal.Length, Petal.Width; OK. Measured: the MCL viewer docks and
  shows "Updating..." for more than two minutes; no column is added and nothing is logged (run p1f, 120 s;
  run p1g fails at the pinned step). With Species and Petal.Length the same dialog finishes in about 2 s
  and adds Cluster (MCL), EmbedX (MCL), EmbedY (MCL), Cluster size (MCL) and Connectivity (MCL). Expected:
  the clustering finishes on any columns the picker offers. Cause seen in the source: since 2808a562b7
  (24 Apr 2026, "ML: Dim Reduction: Skip editor creation for vector metrics") the shared editor builds no
  per-column editors when every picked column is numerical
  (`libraries/ml/src/multi-column-dimensionality-reduction/multi-column-dim-reduction-editor.ts:214`); Reduce
  Dimensionality then passes a vector metric, but `MCLEditor.params` (`libraries/ml/src/MCL/mcl-editor.ts:76`)
  still maps the empty editor list into empty metrics, weights and preprocessing functions, and
  `markovCluster` (`libraries/ml/src/MCL/clustering-view.ts:50`) hands the worker no column at all.
  `clustering.feature`, scenario "MCL on four number columns adds the cluster columns", `@known-failure` on
  the Cluster (MCL) step. The stuck run does not outlive the scenario: closing the MCL viewer, its view or every view
  (`closeAll`, the shell reset after each scenario) terminates its worker within 3 s (probe 30 Oct: live
  MCL workers 1 → 0; `onViewerClosed` in `clustering-view.ts:71`), and the next feature on the same page
  started clean (run rf-mcl2, one worker).
