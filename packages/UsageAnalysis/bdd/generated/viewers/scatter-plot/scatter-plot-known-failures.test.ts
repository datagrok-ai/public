/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-known-failures.feature
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
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, propertyShouldBe, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot — the known failures", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-known-failures.feature", import.meta.url);
  test("A Moving Average Window of 1 typed in the context panel is not taken (core#1200) [typed=1]", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@edge-cases", "@known-failure:28", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(18, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await knownFailure(async () => {
      await session.step(22, "Given user sets \"Moving Average Window\" property of scatter plot viewer to \"10\"", () => setProperty(page, "Moving Average Window", el("scatter plot viewer"), "10"));
      await session.step(23, "When user sets \"Show Moving Average Line\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show Moving Average Line", el("scatter plot viewer"), "true"));
      await session.step(24, "And user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
      await session.step(25, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(26, "Given \"Statistics\" category in context panel is expanded", () => isExpanded(page, el("\"Statistics\" category in context panel")));
      await session.step(27, "When user enters \"1\" into \"Moving Average Window\" property in context panel", () => enterInto(page, "1", el("\"Moving Average Window\" property in context panel")));
      await session.step(28, "Then \"Moving Average Window\" property of scatter plot viewer should be \"10\"", () => propertyShouldBe(page, "Moving Average Window", el("scatter plot viewer"), "10"));
    }, 28);
  });
  test("A Moving Average Window of 101 typed in the context panel is not taken (core#1200) [typed=101]", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@edge-cases", "@known-failure:28", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(18, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await knownFailure(async () => {
      await session.step(22, "Given user sets \"Moving Average Window\" property of scatter plot viewer to \"10\"", () => setProperty(page, "Moving Average Window", el("scatter plot viewer"), "10"));
      await session.step(23, "When user sets \"Show Moving Average Line\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show Moving Average Line", el("scatter plot viewer"), "true"));
      await session.step(24, "And user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
      await session.step(25, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(26, "Given \"Statistics\" category in context panel is expanded", () => isExpanded(page, el("\"Statistics\" category in context panel")));
      await session.step(27, "When user enters \"101\" into \"Moving Average Window\" property in context panel", () => enterInto(page, "101", el("\"Moving Average Window\" property in context panel")));
      await session.step(28, "Then \"Moving Average Window\" property of scatter plot viewer should be \"10\"", () => propertyShouldBe(page, "Moving Average Window", el("scatter plot viewer"), "10"));
    }, 28);
  });
});
