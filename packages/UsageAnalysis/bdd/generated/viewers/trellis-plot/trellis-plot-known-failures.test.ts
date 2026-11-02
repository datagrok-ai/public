/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-known-failures.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot, core#1200]
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
import {addViewerWith, propertyShouldBe, readingIs, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot — the known failures", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-known-failures.feature", import.meta.url);
  test("A Tiles Per Row of 11 typed in the context panel outside its declared bounds is not taken (core#1200)", {tag: ["@viewers", "@realizes:viewers.trellis-plot", "@edge-cases", "@known-failure:29", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","DIS_POP"],["Y Column Names",""],["Viewer Type","Scatter plot"],["Pack Categories","false"]]), [["X Column Names","DIS_POP"],["Y Column Names",""],["Viewer Type","Scatter plot"],["Pack Categories","false"]]);
    await session.step(20, "Then the \"x categories\" reading of trellis plot viewer should be 6", () => readingIs(page, "x categories", el("trellis plot viewer"), 6));
    await knownFailure(async () => {
      await session.step(24, "Given user sets \"Tiles Per Row\" property of trellis plot viewer to \"4\"", () => setProperty(page, "Tiles Per Row", el("trellis plot viewer"), "4"));
      await session.step(25, "When user clicks on settings icon of trellis plot viewer", () => clickOn(page, el("settings icon of trellis plot viewer")));
      await session.step(26, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(27, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
      await session.step(28, "When user enters \"11\" into \"Tiles Per Row\" property in context panel", () => enterInto(page, "11", el("\"Tiles Per Row\" property in context panel")));
      await session.step(29, "Then \"Tiles Per Row\" property of trellis plot viewer should be \"4\"", () => propertyShouldBe(page, "Tiles Per Row", el("trellis plot viewer"), "4"));
    }, 29);
  });
});
