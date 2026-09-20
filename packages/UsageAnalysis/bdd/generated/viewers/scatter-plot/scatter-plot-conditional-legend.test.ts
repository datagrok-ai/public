/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-conditional-legend.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {colorCodedAs, colorConditional, colorOff} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, legendItemColor, legendLists, noErrors, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot legend of a conditionally coloured column", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-conditional-legend.feature", import.meta.url);
  test("A conditionally coloured column set as Color lists one legend entry per rule", {tag: ["@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]));
    await session.step(16, "Then legend of scatter plot viewer should be hidden", () => shouldBe(page, el("legend of scatter plot viewer"), "hidden"));
    await session.step(19, "When user colors \"AGE\" column conditionally:", () => colorConditional(page, "AGE", [["18-45","#00FF00"],["45-89","#FF0000"]]));
    await session.step(22, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
    await session.step(23, "When user sets \"Color\" property of scatter plot viewer to \"AGE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "AGE"));
    await session.step(24, "Then legend of scatter plot viewer should be visible", () => shouldBe(page, el("legend of scatter plot viewer"), "visible"));
    await session.step(25, "And the legend of scatter plot viewer should list 2 items", () => legendLists(page, el("scatter plot viewer"), 2));
    await session.step(26, "And \"18-45\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"18-45\" legend item in legend of scatter plot viewer"), "visible"));
    await session.step(27, "And \"45-89\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"45-89\" legend item in legend of scatter plot viewer"), "visible"));
    await session.step(28, "And the \"18-45\" item in the legend of scatter plot viewer should be colored \"#00FF00\"", () => legendItemColor(page, "18-45", el("scatter plot viewer"), "#00FF00"));
    await session.step(29, "And the \"45-89\" item in the legend of scatter plot viewer should be colored \"#FF0000\"", () => legendItemColor(page, "45-89", el("scatter plot viewer"), "#FF0000"));
    await session.step(30, "When user sets \"Color\" property of scatter plot viewer to \"\"", () => setProperty(page, "Color", el("scatter plot viewer"), ""));
    await session.step(31, "And user removes the coloring of \"AGE\" column", () => colorOff(page, "AGE"));
    await session.step(32, "Then legend of scatter plot viewer should be hidden", () => shouldBe(page, el("legend of scatter plot viewer"), "hidden"));
    await session.step(33, "And no errors should have been logged", () => noErrors(page));
  });
});
