/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tooltips/tooltips-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tooltips]
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
import {check, clickOn, shouldBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {dialogCloses, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hoverArea, noErrors, pickFromContextMenu, pointerAway, setProperties, tooltipColumnsInOrder, tooltipNotColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {checkedInColumnList, toggleInColumnList, uncheckedInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Tooltip edge cases: Edit Tooltip cancelled, every column unticked, a tooltip column removed", () => {
  const session = feature(test, "features/viewers/tooltips/tooltips-edge-cases.feature", import.meta.url);
  test("A change made in Edit Tooltip and cancelled leaves the tooltip as it was", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tooltips"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a box plot viewer", () => addViewer(page, "box plot"));
    await session.step(17, "And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await session.step(19, "When user picks \"Tooltip > Edit...\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("box plot viewer")));
    await session.step(20, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
    await session.step(21, "When user unchecks checkbox in \"Edit Tooltip\" dialog", () => uncheck(page, el("checkbox in \"Edit Tooltip\" dialog")));
    await session.step(22, "And user toggles the \"AGE\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
    await session.step(23, "And user toggles the \"SEX\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "SEX", el("\"Edit Tooltip\" dialog")));
    await session.step(24, "And user toggles the \"WEIGHT\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "WEIGHT", el("\"Edit Tooltip\" dialog")));
    await session.step(25, "And user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
    await session.step(26, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
    await session.step(27, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
    await session.step(28, "Then the tooltip should show columns in the order \"AGE, SEX, WEIGHT\"", () => tooltipColumnsInOrder(page, "AGE, SEX, WEIGHT"));
    await session.step(29, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(32, "When user picks \"Tooltip > Edit...\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("box plot viewer")));
    await session.step(33, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
    await session.step(34, "When user toggles the \"AGE\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
    await session.step(35, "And user toggles the \"RACE\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "RACE", el("\"Edit Tooltip\" dialog")));
    await session.step(36, "Then the \"AGE\" column should not be checked in the column list of \"Edit Tooltip\" dialog", () => uncheckedInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
    await session.step(37, "And the \"RACE\" column should be checked in the column list of \"Edit Tooltip\" dialog", () => checkedInColumnList(page, "RACE", el("\"Edit Tooltip\" dialog")));
    await session.step(38, "When user clicks on CANCEL button in \"Edit Tooltip\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Tooltip\" dialog")));
    await session.step(39, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
    await session.step(40, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
    await session.step(41, "Then the tooltip should show columns in the order \"AGE, SEX, WEIGHT\"", () => tooltipColumnsInOrder(page, "AGE, SEX, WEIGHT"));
    await session.step(42, "And the tooltip should not show columns \"RACE\"", () => tooltipNotColumns(page, "RACE"));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
  });
  test("Every column unticked and confirmed leaves no column in the tooltip", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tooltips"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a box plot viewer", () => addViewer(page, "box plot"));
    await session.step(17, "And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await session.step(19, "When user picks \"Tooltip > Edit...\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("box plot viewer")));
    await session.step(20, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
    await session.step(21, "When user unchecks checkbox in \"Edit Tooltip\" dialog", () => uncheck(page, el("checkbox in \"Edit Tooltip\" dialog")));
    await session.step(22, "And user toggles the \"AGE\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
    await session.step(23, "And user toggles the \"SEX\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "SEX", el("\"Edit Tooltip\" dialog")));
    await session.step(24, "And user toggles the \"WEIGHT\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "WEIGHT", el("\"Edit Tooltip\" dialog")));
    await session.step(25, "And user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
    await session.step(26, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
    await session.step(27, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
    await session.step(28, "Then the tooltip should show columns in the order \"AGE, SEX, WEIGHT\"", () => tooltipColumnsInOrder(page, "AGE, SEX, WEIGHT"));
    await session.step(29, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(46, "When user picks \"Tooltip > Edit...\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("box plot viewer")));
    await session.step(47, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
    await session.step(48, "When user checks checkbox in \"Edit Tooltip\" dialog", () => check(page, el("checkbox in \"Edit Tooltip\" dialog")));
    await session.step(49, "Then the \"RACE\" column should be checked in the column list of \"Edit Tooltip\" dialog", () => checkedInColumnList(page, "RACE", el("\"Edit Tooltip\" dialog")));
    await session.step(50, "When user unchecks checkbox in \"Edit Tooltip\" dialog", () => uncheck(page, el("checkbox in \"Edit Tooltip\" dialog")));
    await session.step(51, "Then the \"AGE\" column should not be checked in the column list of \"Edit Tooltip\" dialog", () => uncheckedInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
    await session.step(52, "And the \"RACE\" column should not be checked in the column list of \"Edit Tooltip\" dialog", () => uncheckedInColumnList(page, "RACE", el("\"Edit Tooltip\" dialog")));
    await session.step(53, "When user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
    await session.step(54, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
    await session.step(55, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
    await session.step(56, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
    await session.step(57, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(58, "And user picks \"Tooltip > Edit...\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("box plot viewer")));
    await session.step(59, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
    await session.step(60, "When user toggles the \"AGE\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
    await session.step(61, "And user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
    await session.step(62, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
    await session.step(63, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
    await session.step(64, "Then the tooltip should show columns in the order \"AGE\"", () => tooltipColumnsInOrder(page, "AGE"));
    await session.step(65, "And no errors should have been logged", () => noErrors(page));
  });
  test("A tooltip column removed from the table is skipped and the others still show", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tooltips"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a box plot viewer", () => addViewer(page, "box plot"));
    await session.step(17, "And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await session.step(19, "When user picks \"Tooltip > Edit...\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("box plot viewer")));
    await session.step(20, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
    await session.step(21, "When user unchecks checkbox in \"Edit Tooltip\" dialog", () => uncheck(page, el("checkbox in \"Edit Tooltip\" dialog")));
    await session.step(22, "And user toggles the \"AGE\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
    await session.step(23, "And user toggles the \"SEX\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "SEX", el("\"Edit Tooltip\" dialog")));
    await session.step(24, "And user toggles the \"WEIGHT\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "WEIGHT", el("\"Edit Tooltip\" dialog")));
    await session.step(25, "And user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
    await session.step(26, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
    await session.step(27, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
    await session.step(28, "Then the tooltip should show columns in the order \"AGE, SEX, WEIGHT\"", () => tooltipColumnsInOrder(page, "AGE, SEX, WEIGHT"));
    await session.step(29, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(68, "When user removes \"SEX\" column", () => removeColumn(page, "SEX"));
    await session.step(69, "And user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
    await session.step(70, "Then the tooltip should show columns in the order \"AGE, WEIGHT\"", () => tooltipColumnsInOrder(page, "AGE, WEIGHT"));
    await session.step(71, "And no errors should have been logged", () => noErrors(page));
  });
});
