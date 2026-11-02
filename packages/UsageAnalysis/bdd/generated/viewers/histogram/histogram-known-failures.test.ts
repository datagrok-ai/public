/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/histogram/histogram-known-failures.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.histogram, core#1200]
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
import {addViewerWith, propertyShouldBe, readingIs, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Histogram known failures", () => {
  const session = feature(test, "features/viewers/histogram/histogram-known-failures.feature", import.meta.url);
  test("A Bins value typed in the context panel outside its declared bounds is not taken (core#1200)", {tag: ["@viewers", "@realizes:viewers.histogram", "@edge-cases", "@known-failure:23", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"]]), [["Value","AGE"]]);
    await session.step(14, "Then histogram viewer should show 1000 rows", () => showsRows(page, el("histogram viewer"), 1000));
    await session.step(15, "And the \"bins shown\" reading of histogram viewer should be 20", () => readingIs(page, "bins shown", el("histogram viewer"), 20));
    await knownFailure(async () => {
      await session.step(19, "When user clicks on settings icon of histogram viewer", () => clickOn(page, el("settings icon of histogram viewer")));
      await session.step(20, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(21, "Given \"X Axis\" category in context panel is expanded", () => isExpanded(page, el("\"X Axis\" category in context panel")));
      await session.step(22, "When user enters \"201\" into \"Bins\" property in context panel", () => enterInto(page, "201", el("\"Bins\" property in context panel")));
      await session.step(23, "Then \"Bins\" property of histogram viewer should be \"20\"", () => propertyShouldBe(page, "Bins", el("histogram viewer"), "20"));
    }, 23);
  });
});
