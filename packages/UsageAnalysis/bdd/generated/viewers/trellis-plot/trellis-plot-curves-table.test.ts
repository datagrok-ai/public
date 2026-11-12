/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-curves-table.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
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
import {clickOn, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset, packageInstalled} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, boundTable, noErrors, readingIs, readingReads, reportsNoError, setProperties} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInnerViewer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot on another table, with the Multi Curve viewer inside", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-curves-table.feature", import.meta.url);
  test("Trellis plot on another table, with the Multi Curve viewer inside", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(17, "And user opens curves dataset", () => openDataset(page, ds("curves")));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(23, "Then the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
    await session.step(24, "And trellis plot viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("trellis plot viewer"), "demog-1000"));
    await run.scenario("Set to the curves table in the context panel, the trellis works on it and takes the Multi Curve viewer", async () => {
      await session.step(27, "When user clicks on settings icon of trellis plot viewer", () => clickOn(page, el("settings icon of trellis plot viewer")));
      await session.step(28, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(29, "When user selects \"curves\" in \"Table\" property in context panel", () => selectIn(page, "curves", el("\"Table\" property in context panel")));
      await session.step(30, "Then trellis plot viewer should be bound to table \"curves\"", () => boundTable(page, el("trellis plot viewer"), "curves"));
      await session.step(31, "And the \"table rows\" reading of trellis plot viewer should be 30", () => readingIs(page, "table rows", el("trellis plot viewer"), 30));
      await session.step(32, "When user picks \"Curves\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Curves", el("trellis plot viewer")));
      await session.step(33, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"MultiCurveViewer\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "MultiCurveViewer"));
      await session.step(34, "And the \"cells drawn\" reading of trellis plot viewer should be 1", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 1));
      await session.step(35, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
      await session.step(36, "And trellis plot viewer should report no error", () => reportsNoError(page, el("trellis plot viewer")));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Set back to demog-1000, the canonical grid comes back", async () => {
      await session.step(40, "When user selects \"demog-1000\" in \"Table\" property in context panel", () => selectIn(page, "demog-1000", el("\"Table\" property in context panel")));
      await session.step(41, "And user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
      await session.step(45, "Then trellis plot viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("trellis plot viewer"), "demog-1000"));
      await session.step(46, "And the \"table rows\" reading of trellis plot viewer should be 1000", () => readingIs(page, "table rows", el("trellis plot viewer"), 1000));
      await session.step(47, "And the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(48, "And the \"rows shown\" reading of trellis plot viewer should be 1000", () => readingIs(page, "rows shown", el("trellis plot viewer"), 1000));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
