/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/grid-reorder-and-hide-columns.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dragAreaToArea, hasNoArea, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingExcludes, readingIncludes, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Reorder columns by dragging and hide the ones you do not need", () => {
  const session = feature(test, "features/guides/grid-reorder-and-hide-columns.feature", import.meta.url);
  test("Drag a column header to a new place, then hide two columns", {tag: ["@guide", "@help:visualize/viewers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(13, "When user drags the \"header HEIGHT\" area of grid to the \"header AGE\" area", () => dragAreaToArea(page, "header HEIGHT", el("grid"), "header AGE"));
    await session.step(14, "Then the \"column order\" reading of grid should include the text \"AGE, HEIGHT, SEX\"", () => readingIncludes(page, "column order", el("grid"), "AGE, HEIGHT, SEX"));
    await session.step(15, "When user picks \"Hide\" from the context menu of the \"header WEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Hide", "header WEIGHT", el("grid")));
    await session.step(16, "Then grid should not have a \"header WEIGHT\" area", () => hasNoArea(page, el("grid"), "header WEIGHT"));
    await session.step(17, "When user picks \"Order or Hide Columns...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Order or Hide Columns...", "header AGE", el("grid")));
    await session.step(18, "Then Order or Hide Columns dialog should be visible", () => shouldBe(page, el("Order or Hide Columns dialog"), "visible"));
    await session.step(19, "When user toggles the \"STARTED\" column in the column list of Order or Hide Columns dialog", () => toggleInColumnList(page, "STARTED", el("Order or Hide Columns dialog")));
    await session.step(20, "And user clicks on CLOSE button in Order or Hide Columns dialog", () => clickOn(page, el("CLOSE button in Order or Hide Columns dialog")));
    await session.step(21, "Then grid should not have a \"header STARTED\" area", () => hasNoArea(page, el("grid"), "header STARTED"));
    await session.step(22, "And the \"column order\" reading of grid should not include the text \"WEIGHT\"", () => readingExcludes(page, "column order", el("grid"), "WEIGHT"));
  });
});
