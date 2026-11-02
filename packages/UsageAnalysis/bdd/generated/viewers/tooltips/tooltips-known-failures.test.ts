/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tooltips/tooltips-known-failures.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tooltips, core#1215]
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
import {check, clickOn, enterInto, shouldBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {dialogCloses, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hoverArea, pickFromAreaContextMenu, pickFromContextMenu, pointerAway, setProperties, tooltipColumnsInOrder} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {scrollGridTo, toggleInColumnList, uncheckedInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Tooltips — the known failures", () => {
  const session = feature(test, "features/viewers/tooltips/tooltips-known-failures.feature", import.meta.url);
  test("A tooltip column renamed is shown under its new name", {tag: ["@viewers", "@realizes:viewers.tooltips", "@edge-cases", "@known-failure:42", "@realizes:core#1215"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user sets properties of grid:", () => setProperties(page, el("grid"), [["Show Tooltip","inherit from table"],["Show Column Names","Always"],["Show Visible Columns In Tooltip","true"]]), [["Show Tooltip","inherit from table"],["Show Column Names","Always"],["Show Visible Columns In Tooltip","true"]]);
    await session.step(19, "And user adds a box plot viewer", () => addViewer(page, "box plot"));
    await session.step(20, "And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await knownFailure(async () => {
      await session.step(25, "When user picks \"Tooltip > Edit...\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("box plot viewer")));
      await session.step(26, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
      await session.step(27, "When user checks checkbox in \"Edit Tooltip\" dialog", () => check(page, el("checkbox in \"Edit Tooltip\" dialog")));
      await session.step(28, "And user unchecks checkbox in \"Edit Tooltip\" dialog", () => uncheck(page, el("checkbox in \"Edit Tooltip\" dialog")));
      await session.step(29, "Then the \"WEIGHT\" column should not be checked in the column list of \"Edit Tooltip\" dialog", () => uncheckedInColumnList(page, "WEIGHT", el("\"Edit Tooltip\" dialog")));
      await session.step(30, "When user toggles the \"AGE\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
      await session.step(31, "And user toggles the \"WEIGHT\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "WEIGHT", el("\"Edit Tooltip\" dialog")));
      await session.step(32, "And user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
      await session.step(33, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
      await session.step(34, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(35, "Then the tooltip should show columns in the order \"AGE, WEIGHT\"", () => tooltipColumnsInOrder(page, "AGE, WEIGHT"));
      await session.step(36, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(37, "And user scrolls the grid to the \"WEIGHT\" column", () => scrollGridTo(page, "WEIGHT"));
      await session.step(38, "And user picks \"Column Properties...\" from the context menu of the \"header WEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header WEIGHT", el("grid")));
      await session.step(39, "And user enters \"MASS\" into \"New name:\" input in \"WEIGHT\" dialog", () => enterInto(page, "MASS", el("\"New name:\" input in \"WEIGHT\" dialog")));
      await session.step(40, "And user clicks on OK button in \"WEIGHT\" dialog", () => clickOn(page, el("OK button in \"WEIGHT\" dialog")));
      await session.step(41, "And user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(42, "Then the tooltip should show columns in the order \"AGE, MASS\"", () => tooltipColumnsInOrder(page, "AGE, MASS"));
    }, 42);
  });
});
