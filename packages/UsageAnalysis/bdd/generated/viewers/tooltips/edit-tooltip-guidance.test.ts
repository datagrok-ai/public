/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tooltips/edit-tooltip-guidance.feature
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
import {check, clickOn, hoverOver, selectIn, shouldBe, shouldContainText, shouldHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {dialogCloses, openDataset, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, closeContextMenu, hoverArea, menuListsDisabled, menuListsEnabled, noErrors, openContextMenu, pickFromContextMenu, pointerAway, propertyShouldBe, setProperties, tooltipColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {checkedInColumnList, toggleInColumnList, uncheckedInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Edit Tooltip dialog says what it edits and what the tooltip will show", () => {
  const session = feature(test, "features/viewers/tooltips/edit-tooltip-guidance.feature", import.meta.url);
  test("The Edit Tooltip dialog says what it edits and what the tooltip will show", {tag: ["@journey", "@viewers", "@realizes:viewers.tooltips"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(23, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(24, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await run.scenario("Use tooltip and Show column names explain themselves on hover", async () => {
      await session.step(28, "When user picks \"Tooltip > Edit...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("scatter plot viewer")));
      await session.step(29, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
      await session.step(30, "And \"Use tooltip\" input in \"Edit Tooltip\" dialog should have value \"Table (shared)\"", () => shouldHaveValue(page, el("\"Use tooltip\" input in \"Edit Tooltip\" dialog"), "Table (shared)"));
      await session.step(31, "When user hovers over \"Use tooltip\" label in \"Edit Tooltip\" dialog", () => hoverOver(page, el("\"Use tooltip\" label in \"Edit Tooltip\" dialog")));
      await session.step(32, "Then tooltip should contain text \"every viewer whose Show Tooltip is\"", () => shouldContainText(page, el("tooltip"), "every viewer whose Show Tooltip is"));
      await session.step(33, "And tooltip should contain text \"Scatter plot only\"", () => shouldContainText(page, el("tooltip"), "Scatter plot only"));
      await session.step(34, "When user hovers over \"Show column names\" label in \"Edit Tooltip\" dialog", () => hoverOver(page, el("\"Show column names\" label in \"Edit Tooltip\" dialog")));
      await session.step(35, "Then tooltip should contain text \"only where values could be mixed up\"", () => shouldContainText(page, el("tooltip"), "only where values could be mixed up"));
      await session.step(36, "When user clicks on CANCEL button in \"Edit Tooltip\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Tooltip\" dialog")));
      await session.step(37, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The status line counts the checked columns, warns past the table's 7 and on none, and each tooltip keeps its checks", async () => {
      await session.step(41, "When user picks \"Tooltip > Edit...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("scatter plot viewer")));
      await session.step(42, "Then \"Edit Tooltip\" dialog should contain text \"Shows the first 7 of 11 checked\"", () => shouldContainText(page, el("\"Edit Tooltip\" dialog"), "Shows the first 7 of 11 checked"));
      await session.step(43, "When user unchecks checkbox in \"Edit Tooltip\" dialog", () => uncheck(page, el("checkbox in \"Edit Tooltip\" dialog")));
      await session.step(44, "Then \"Edit Tooltip\" dialog should contain text \"Nothing checked: no row tooltip\"", () => shouldContainText(page, el("\"Edit Tooltip\" dialog"), "Nothing checked: no row tooltip"));
      await session.step(45, "When user toggles the \"AGE\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
      await session.step(46, "Then \"Edit Tooltip\" dialog should contain text \"1 of 11 checked\"", () => shouldContainText(page, el("\"Edit Tooltip\" dialog"), "1 of 11 checked"));
      await session.step(47, "When user selects \"Scatter plot only\" in \"Use tooltip\" input in \"Edit Tooltip\" dialog", () => selectIn(page, "Scatter plot only", el("\"Use tooltip\" input in \"Edit Tooltip\" dialog")));
      await session.step(48, "Then \"Edit Tooltip\" dialog should contain text \"Nothing checked: only X and Y show\"", () => shouldContainText(page, el("\"Edit Tooltip\" dialog"), "Nothing checked: only X and Y show"));
      await session.step(49, "When user checks checkbox in \"Edit Tooltip\" dialog", () => check(page, el("checkbox in \"Edit Tooltip\" dialog")));
      await session.step(50, "Then \"Edit Tooltip\" dialog should contain text \"11 of 11 checked\"", () => shouldContainText(page, el("\"Edit Tooltip\" dialog"), "11 of 11 checked"));
      await session.step(51, "When user selects \"Table (shared)\" in \"Use tooltip\" input in \"Edit Tooltip\" dialog", () => selectIn(page, "Table (shared)", el("\"Use tooltip\" input in \"Edit Tooltip\" dialog")));
      await session.step(52, "Then \"Edit Tooltip\" dialog should contain text \"1 of 11 checked\"", () => shouldContainText(page, el("\"Edit Tooltip\" dialog"), "1 of 11 checked"));
      await session.step(53, "And the \"AGE\" column should be checked in the column list of \"Edit Tooltip\" dialog", () => checkedInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
      await session.step(54, "And the \"SEX\" column should not be checked in the column list of \"Edit Tooltip\" dialog", () => uncheckedInColumnList(page, "SEX", el("\"Edit Tooltip\" dialog")));
      await session.step(55, "When user clicks on CANCEL button in \"Edit Tooltip\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Tooltip\" dialog")));
      await session.step(56, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("\"Scatter plot only\" with every column checked shows all 11, where the table's tooltip shows 7", async () => {
      await session.step(60, "When user picks \"Tooltip > Edit...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("scatter plot viewer")));
      await session.step(61, "And user selects \"Scatter plot only\" in \"Use tooltip\" input in \"Edit Tooltip\" dialog", () => selectIn(page, "Scatter plot only", el("\"Use tooltip\" input in \"Edit Tooltip\" dialog")));
      await session.step(62, "And user checks checkbox in \"Edit Tooltip\" dialog", () => check(page, el("checkbox in \"Edit Tooltip\" dialog")));
      await session.step(63, "And user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
      await session.step(64, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
      await session.step(65, "And \"Show Tooltip\" property of scatter plot viewer should be \"show custom tooltip\"", () => propertyShouldBe(page, "Show Tooltip", el("scatter plot viewer"), "show custom tooltip"));
      await session.step(66, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(67, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(68, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(69, "And user picks \"Tooltip > Edit...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("scatter plot viewer")));
      await session.step(70, "And user selects \"Table (shared)\" in \"Use tooltip\" input in \"Edit Tooltip\" dialog", () => selectIn(page, "Table (shared)", el("\"Use tooltip\" input in \"Edit Tooltip\" dialog")));
      await session.step(71, "Then \"Edit Tooltip\" dialog should contain text \"Shows the first 7 of 11 checked\"", () => shouldContainText(page, el("\"Edit Tooltip\" dialog"), "Shows the first 7 of 11 checked"));
      await session.step(72, "When user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
      await session.step(73, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
      await session.step(74, "And \"Show Tooltip\" property of scatter plot viewer should be \"inherit from table\"", () => propertyShouldBe(page, "Show Tooltip", el("scatter plot viewer"), "inherit from table"));
      await session.step(75, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(76, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(77, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(78, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Design custom tooltip... designs the scatter plot's own tooltip once Use tooltip says so, and the dialog then says OK replaces it", async () => {
      await session.step(81, "When user picks \"Tooltip > Edit...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("scatter plot viewer")));
      await session.step(82, "And user selects \"Scatter plot only\" in \"Use tooltip\" input in \"Edit Tooltip\" dialog", () => selectIn(page, "Scatter plot only", el("\"Use tooltip\" input in \"Edit Tooltip\" dialog")));
      await session.step(83, "And user clicks on \"Design custom tooltip...\" label in \"Edit Tooltip\" dialog", () => clickOn(page, el("\"Design custom tooltip...\" label in \"Edit Tooltip\" dialog")));
      await session.step(84, "And user clicks on \"CLOSE AND APPLY\" button", () => clickOn(page, el("\"CLOSE AND APPLY\" button")));
      await session.step(85, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
      await session.step(86, "And \"Show Tooltip\" property of scatter plot viewer should be \"show custom tooltip\"", () => propertyShouldBe(page, "Show Tooltip", el("scatter plot viewer"), "show custom tooltip"));
      await session.step(87, "When user picks \"Tooltip > Edit...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("scatter plot viewer")));
      await session.step(88, "Then \"Use tooltip\" input in \"Edit Tooltip\" dialog should have value \"Scatter plot only\"", () => shouldHaveValue(page, el("\"Use tooltip\" input in \"Edit Tooltip\" dialog"), "Scatter plot only"));
      await session.step(89, "And \"Edit Tooltip\" dialog should contain text \"A designed tooltip is in use: OK replaces it\"", () => shouldContainText(page, el("\"Edit Tooltip\" dialog"), "A designed tooltip is in use: OK replaces it"));
      await session.step(90, "When user clicks on CANCEL button in \"Edit Tooltip\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Tooltip\" dialog")));
      await session.step(91, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Remove Group Tooltip is disabled until a viewer is the group tooltip, which Reset group tooltip removes on OK", async () => {
      await session.step(95, "When user opens the context menu of scatter plot viewer", () => openContextMenu(page, el("scatter plot viewer")));
      await session.step(96, "Then the open menu should list \"Tooltip > Remove Group Tooltip\" as disabled", () => menuListsDisabled(page, "Tooltip > Remove Group Tooltip"));
      await session.step(97, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(98, "And user picks \"Tooltip > Use as Group Tooltip\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Use as Group Tooltip", el("scatter plot viewer")));
      await session.step(99, "And user opens the context menu of scatter plot viewer", () => openContextMenu(page, el("scatter plot viewer")));
      await session.step(100, "Then the open menu should list \"Tooltip > Remove Group Tooltip\" as enabled", () => menuListsEnabled(page, "Tooltip > Remove Group Tooltip"));
      await session.step(101, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(102, "And user picks \"Tooltip > Edit...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("scatter plot viewer")));
      await session.step(103, "And user clicks on \"Reset group tooltip\" label in \"Edit Tooltip\" dialog", () => clickOn(page, el("\"Reset group tooltip\" label in \"Edit Tooltip\" dialog")));
      await session.step(104, "Then \"Edit Tooltip\" dialog should contain text \"The group tooltip is removed on OK\"", () => shouldContainText(page, el("\"Edit Tooltip\" dialog"), "The group tooltip is removed on OK"));
      await session.step(105, "When user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
      await session.step(106, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
      await session.step(107, "When user opens the context menu of scatter plot viewer", () => openContextMenu(page, el("scatter plot viewer")));
      await session.step(108, "Then the open menu should list \"Tooltip > Remove Group Tooltip\" as disabled", () => menuListsDisabled(page, "Tooltip > Remove Group Tooltip"));
      await session.step(109, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(110, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
