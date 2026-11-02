/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tree/tree-on-click-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.tree, public#3249, public#3298]
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
import {clickOn, isExpanded, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, filterIsExactlyCategory, filterPasses, noneOfFiltered, noneOfSelected, noneSelected, resetFilter, selectWhereIs, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, clickAreaHolding, hasArea, lessHighlight, moreHighlight, noErrors, propertyShouldBe, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Tree clicks: On Click = Filter and Select, Shift and Control, and the grid's selection drawn on the tree", () => {
  const session = feature(test, "features/viewers/tree/tree-on-click-filter.feature", import.meta.url);
  test("The grid's selection is drawn on the tree; a click filters to the branch, Shift or Control adds a branch, a plain click starts over, and Select selects again (public#3249, public#3298)", {tag: ["@viewers", "@realizes:charts.viewer.tree", "@realizes:public#3249", "@realizes:public#3298"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(21, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(22, "When user clicks on first \"Tree\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Tree\" card in \"Add Viewer\" dialog")));
    await session.step(23, "Then tree viewer should be visible", () => shouldBe(page, el("tree viewer"), "visible"));
    await session.step(24, "When user clicks on grid", () => clickOn(page, el("grid")));
    await session.step(25, "And user clicks on settings icon of tree viewer", () => clickOn(page, el("settings icon of tree viewer")));
    await session.step(26, "And user clicks on \"...\" button in \"Hierarchy\" property in context panel", () => clickOn(page, el("\"...\" button in \"Hierarchy\" property in context panel")));
    await session.step(27, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(28, "When user clicks on \"None\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"None\" link in \"Select columns...\" dialog")));
    await session.step(29, "Then \"0 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"0 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(30, "When user toggles the \"CONTROL\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "CONTROL", el("\"Select columns...\" dialog")));
    await session.step(31, "And user toggles the \"SEX\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "SEX", el("\"Select columns...\" dialog")));
    await session.step(32, "And user toggles the \"RACE\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "RACE", el("\"Select columns...\" dialog")));
    await session.step(33, "Then \"3 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"3 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(34, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(35, "Then \"Select columns...\" dialog should be absent", () => shouldBe(page, el("\"Select columns...\" dialog"), "absent"));
    await session.step(36, "And \"Hierarchy\" property of tree viewer should be \"CONTROL, SEX, RACE\"", () => propertyShouldBe(page, "Hierarchy", el("tree viewer"), "CONTROL, SEX, RACE"));
    await session.step(37, "And tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(41, "When user takes a snapshot of tree viewer", () => takeSnapshot(page, el("tree viewer")));
    await session.step(42, "And user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
    await session.step(43, "Then 15 rows should be selected", () => selectedRowCount(page, 15));
    await session.step(44, "And tree viewer should show more selection highlight than before", () => moreHighlight(page, el("tree viewer")));
    await session.step(45, "When user takes a snapshot of tree viewer", () => takeSnapshot(page, el("tree viewer")));
    await session.step(46, "And user clears the row selection", () => clearSelection(page));
    await session.step(47, "Then no rows should be selected", () => noneSelected(page));
    await session.step(48, "And tree viewer should show less selection highlight than before", () => lessHighlight(page, el("tree viewer")));
    await session.step(49, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
    await session.step(50, "When user selects \"Filter\" in \"On Click\" property in context panel", () => selectIn(page, "Filter", el("\"On Click\" property in context panel")));
    await session.step(51, "Then \"On Click\" property of tree viewer should be \"Filter\"", () => propertyShouldBe(page, "On Click", el("tree viewer"), "Filter"));
    await session.step(52, "And 1000 rows should pass the filter", () => filterPasses(page, 1000));
    await session.step(53, "When user clicks on the \"branch All | false | F | Asian\" area of tree viewer", () => clickArea(page, "branch All | false | F | Asian", el("tree viewer")));
    await session.step(54, "Then 7 rows should pass the filter", () => filterPasses(page, 7));
    await session.step(55, "And no rows where \"SEX\" is \"M\" should pass the filter", () => noneOfFiltered(page, "SEX", "M"));
    await session.step(56, "And no rows where \"RACE\" is \"Caucasian\" should pass the filter", () => noneOfFiltered(page, "RACE", "Caucasian"));
    await session.step(57, "And no rows should be selected", () => noneSelected(page));
    await session.step(58, "When user clicks on the \"branch All | false | M | Asian\" area of tree viewer holding Shift", () => clickAreaHolding(page, "branch All | false | M | Asian", el("tree viewer"), "Shift"));
    await session.step(59, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
    await session.step(60, "And the filter should pass exactly the rows where \"RACE\" is \"Asian\"", () => filterIsExactlyCategory(page, "RACE", "Asian"));
    await session.step(61, "When user clicks on the \"branch All | true | F | Caucasian\" area of tree viewer", () => clickArea(page, "branch All | true | F | Caucasian", el("tree viewer")));
    await session.step(62, "Then 5 rows should pass the filter", () => filterPasses(page, 5));
    await session.step(63, "And no rows where \"CONTROL\" is \"false\" should pass the filter", () => noneOfFiltered(page, "CONTROL", "false"));
    await session.step(64, "And no rows where \"RACE\" is \"Asian\" should pass the filter", () => noneOfFiltered(page, "RACE", "Asian"));
    await session.step(65, "And no rows should be selected", () => noneSelected(page));
    await session.step(66, "When user clicks on the \"branch All | false | F | Asian\" area of tree viewer", () => clickArea(page, "branch All | false | F | Asian", el("tree viewer")));
    await session.step(67, "Then 7 rows should pass the filter", () => filterPasses(page, 7));
    await session.step(68, "When user clicks on the \"branch All | false | M | Asian\" area of tree viewer holding Control", () => clickAreaHolding(page, "branch All | false | M | Asian", el("tree viewer"), "Control"));
    await session.step(69, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
    await session.step(70, "And the filter should pass exactly the rows where \"RACE\" is \"Asian\"", () => filterIsExactlyCategory(page, "RACE", "Asian"));
    await session.step(71, "And no rows should be selected", () => noneSelected(page));
    await session.step(72, "When user selects \"Select\" in \"On Click\" property in context panel", () => selectIn(page, "Select", el("\"On Click\" property in context panel")));
    await session.step(73, "Then \"On Click\" property of tree viewer should be \"Select\"", () => propertyShouldBe(page, "On Click", el("tree viewer"), "Select"));
    await session.step(74, "When user resets the filter", () => resetFilter(page));
    await session.step(75, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
    await session.step(76, "When user clicks on the \"branch All | false | F | Asian\" area of tree viewer", () => clickArea(page, "branch All | false | F | Asian", el("tree viewer")));
    await session.step(77, "Then 7 rows should be selected", () => selectedRowCount(page, 7));
    await session.step(78, "And no rows where \"SEX\" is \"M\" should be selected", () => noneOfSelected(page, "SEX", "M"));
    await session.step(79, "When user clicks on the \"branch All | false | M | Asian\" area of tree viewer holding Control", () => clickAreaHolding(page, "branch All | false | M | Asian", el("tree viewer"), "Control"));
    await session.step(80, "Then 15 rows should be selected", () => selectedRowCount(page, 15));
    await session.step(81, "And no errors should have been logged", () => noErrors(page));
  });
  test("A first click made with Shift filters to the branch as a plain click does", {tag: ["@viewers", "@realizes:charts.viewer.tree"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(21, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(22, "When user clicks on first \"Tree\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Tree\" card in \"Add Viewer\" dialog")));
    await session.step(23, "Then tree viewer should be visible", () => shouldBe(page, el("tree viewer"), "visible"));
    await session.step(24, "When user clicks on grid", () => clickOn(page, el("grid")));
    await session.step(25, "And user clicks on settings icon of tree viewer", () => clickOn(page, el("settings icon of tree viewer")));
    await session.step(26, "And user clicks on \"...\" button in \"Hierarchy\" property in context panel", () => clickOn(page, el("\"...\" button in \"Hierarchy\" property in context panel")));
    await session.step(27, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(28, "When user clicks on \"None\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"None\" link in \"Select columns...\" dialog")));
    await session.step(29, "Then \"0 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"0 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(30, "When user toggles the \"CONTROL\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "CONTROL", el("\"Select columns...\" dialog")));
    await session.step(31, "And user toggles the \"SEX\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "SEX", el("\"Select columns...\" dialog")));
    await session.step(32, "And user toggles the \"RACE\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "RACE", el("\"Select columns...\" dialog")));
    await session.step(33, "Then \"3 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"3 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(34, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(35, "Then \"Select columns...\" dialog should be absent", () => shouldBe(page, el("\"Select columns...\" dialog"), "absent"));
    await session.step(36, "And \"Hierarchy\" property of tree viewer should be \"CONTROL, SEX, RACE\"", () => propertyShouldBe(page, "Hierarchy", el("tree viewer"), "CONTROL, SEX, RACE"));
    await session.step(37, "And tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(84, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
    await session.step(85, "When user selects \"Filter\" in \"On Click\" property in context panel", () => selectIn(page, "Filter", el("\"On Click\" property in context panel")));
    await session.step(86, "Then \"On Click\" property of tree viewer should be \"Filter\"", () => propertyShouldBe(page, "On Click", el("tree viewer"), "Filter"));
    await session.step(87, "And 1000 rows should pass the filter", () => filterPasses(page, 1000));
    await session.step(88, "When user clicks on the \"branch All | false | F | Asian\" area of tree viewer holding Shift", () => clickAreaHolding(page, "branch All | false | F | Asian", el("tree viewer"), "Shift"));
    await session.step(89, "Then 7 rows should pass the filter", () => filterPasses(page, 7));
    await session.step(90, "And no errors should have been logged", () => noErrors(page));
  });
});
