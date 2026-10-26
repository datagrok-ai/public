/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/correlation-plot/correlation-plot-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.correlation-plot]
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
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hasArea, hasNoArea, noErrors, pickFromAreaContextMenu, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Correlation plot regressions of closed issues", () => {
  const session = feature(test, "features/viewers/correlation-plot/correlation-plot-regressions.feature", import.meta.url);
  test("Hiding a column in Order or Hide Columns keeps the row names (public#3492)", {tag: ["@viewers", "@realizes:viewers.correlation-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a correlation plot viewer", () => addViewer(page, "correlation plot"));
    await session.step(12, "Then the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
    await session.step(13, "And correlation plot viewer should have a \"header WEIGHT\" area", () => hasArea(page, el("correlation plot viewer"), "header WEIGHT"));
    await session.step(14, "And correlation plot viewer should have a \"row header HEIGHT\" area", () => hasArea(page, el("correlation plot viewer"), "row header HEIGHT"));
    await session.step(15, "When user picks \"Grid > Order or Hide Columns...\" from the context menu of the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => pickFromAreaContextMenu(page, "Grid > Order or Hide Columns...", "cell HEIGHT x AGE", el("correlation plot viewer")));
    await session.step(16, "Then Order or Hide Columns dialog should be visible", () => shouldBe(page, el("Order or Hide Columns dialog"), "visible"));
    await session.step(17, "When user toggles the \"WEIGHT\" column in the column list of Order or Hide Columns dialog", () => toggleInColumnList(page, "WEIGHT", el("Order or Hide Columns dialog")));
    await session.step(18, "And user clicks on CLOSE button in Order or Hide Columns dialog", () => clickOn(page, el("CLOSE button in Order or Hide Columns dialog")));
    await session.step(19, "Then Order or Hide Columns dialog should be hidden", () => shouldBe(page, el("Order or Hide Columns dialog"), "hidden"));
    await session.step(20, "And correlation plot viewer should not have a \"header WEIGHT\" area", () => hasNoArea(page, el("correlation plot viewer"), "header WEIGHT"));
    await session.step(21, "And correlation plot viewer should have a \"pinned band\" area", () => hasArea(page, el("correlation plot viewer"), "pinned band"));
    await session.step(22, "And correlation plot viewer should have a \"row header HEIGHT\" area", () => hasArea(page, el("correlation plot viewer"), "row header HEIGHT"));
    await session.step(23, "And correlation plot viewer should have a \"row header AGE\" area", () => hasArea(page, el("correlation plot viewer"), "row header AGE"));
    await session.step(24, "And correlation plot viewer should have a \"row header WEIGHT\" area", () => hasArea(page, el("correlation plot viewer"), "row header WEIGHT"));
    await session.step(25, "When user picks \"Grid > Order or Hide Columns...\" from the context menu of the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => pickFromAreaContextMenu(page, "Grid > Order or Hide Columns...", "cell HEIGHT x AGE", el("correlation plot viewer")));
    await session.step(26, "And user toggles the \"WEIGHT\" column in the column list of Order or Hide Columns dialog", () => toggleInColumnList(page, "WEIGHT", el("Order or Hide Columns dialog")));
    await session.step(27, "And user clicks on CLOSE button in Order or Hide Columns dialog", () => clickOn(page, el("CLOSE button in Order or Hide Columns dialog")));
    await session.step(28, "Then correlation plot viewer should have a \"header WEIGHT\" area", () => hasArea(page, el("correlation plot viewer"), "header WEIGHT"));
    await session.step(29, "And correlation plot viewer should have a \"row header HEIGHT\" area", () => hasArea(page, el("correlation plot viewer"), "row header HEIGHT"));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
  });
});
