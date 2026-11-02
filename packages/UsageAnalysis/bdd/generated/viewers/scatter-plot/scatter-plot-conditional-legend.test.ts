/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-conditional-legend.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, core#1200]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {colorCodedAs, colorConditional, colorOff} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, legendItemColor, legendLists, noErrors, propertyShouldBe, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot legend of a conditionally coloured column", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-conditional-legend.feature", import.meta.url);
  test("A conditionally coloured column set as Color lists one legend entry per rule", {tag: ["@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(21, "Then legend of scatter plot viewer should be hidden", () => shouldBe(page, el("legend of scatter plot viewer"), "hidden"));
    await session.step(24, "When user colors \"AGE\" column conditionally:", () => colorConditional(page, "AGE", [["18-45","#00FF00"],["45-89","#FF0000"]]), [["18-45","#00FF00"],["45-89","#FF0000"]]);
    await session.step(27, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
    await session.step(28, "When user sets \"Color\" property of scatter plot viewer to \"AGE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "AGE"));
    await session.step(29, "Then legend of scatter plot viewer should be visible", () => shouldBe(page, el("legend of scatter plot viewer"), "visible"));
    await session.step(30, "And the legend of scatter plot viewer should list 2 items", () => legendLists(page, el("scatter plot viewer"), 2));
    await session.step(31, "And \"18-45\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"18-45\" legend item in legend of scatter plot viewer"), "visible"));
    await session.step(32, "And \"45-89\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"45-89\" legend item in legend of scatter plot viewer"), "visible"));
    await session.step(33, "And the \"18-45\" item in the legend of scatter plot viewer should be colored \"#00FF00\"", () => legendItemColor(page, "18-45", el("scatter plot viewer"), "#00FF00"));
    await session.step(34, "And the \"45-89\" item in the legend of scatter plot viewer should be colored \"#FF0000\"", () => legendItemColor(page, "45-89", el("scatter plot viewer"), "#FF0000"));
    await session.step(35, "When user sets \"Color\" property of scatter plot viewer to \"\"", () => setProperty(page, "Color", el("scatter plot viewer"), ""));
    await session.step(36, "And user removes the coloring of \"AGE\" column", () => colorOff(page, "AGE"));
    await session.step(37, "Then legend of scatter plot viewer should be hidden", () => shouldBe(page, el("legend of scatter plot viewer"), "hidden"));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
  test("A Moving Average Window of 1 typed in the context panel is not taken (core#1200) [typed=1]", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:core#1200", "@known-failure:47"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(21, "Then legend of scatter plot viewer should be hidden", () => shouldBe(page, el("legend of scatter plot viewer"), "hidden"));
    await knownFailure(async () => {
      await session.step(42, "When user sets \"Show Moving Average Line\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show Moving Average Line", el("scatter plot viewer"), "true"));
      await session.step(43, "And user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
      await session.step(44, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(45, "Given \"Statistics\" category in context panel is expanded", () => isExpanded(page, el("\"Statistics\" category in context panel")));
      await session.step(46, "When user enters \"1\" into \"Moving Average Window\" property in context panel", () => enterInto(page, "1", el("\"Moving Average Window\" property in context panel")));
      await session.step(47, "Then \"Moving Average Window\" property of scatter plot viewer should be \"10\"", () => propertyShouldBe(page, "Moving Average Window", el("scatter plot viewer"), "10"));
    }, 47);
  });
  test("A Moving Average Window of 101 typed in the context panel is not taken (core#1200) [typed=101]", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:core#1200", "@known-failure:47"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(21, "Then legend of scatter plot viewer should be hidden", () => shouldBe(page, el("legend of scatter plot viewer"), "hidden"));
    await knownFailure(async () => {
      await session.step(42, "When user sets \"Show Moving Average Line\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show Moving Average Line", el("scatter plot viewer"), "true"));
      await session.step(43, "And user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
      await session.step(44, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(45, "Given \"Statistics\" category in context panel is expanded", () => isExpanded(page, el("\"Statistics\" category in context panel")));
      await session.step(46, "When user enters \"101\" into \"Moving Average Window\" property in context panel", () => enterInto(page, "101", el("\"Moving Average Window\" property in context panel")));
      await session.step(47, "Then \"Moving Average Window\" property of scatter plot viewer should be \"10\"", () => propertyShouldBe(page, "Moving Average Window", el("scatter plot viewer"), "10"));
    }, 47);
  });
});
