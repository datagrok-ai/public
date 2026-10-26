/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tree/tree-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.tree]
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
import {clearSelection, filterIsExactlyCategory, filterPasses, noneOfSelected, noneSelected, resetFilter, selectWhereIs, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, clickAreaHolding, hasArea, lessHighlight, moreHighlight, noErrors, propertyShouldBe, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Tree regressions of closed issues", () => {
  const session = feature(test, "features/viewers/tree/tree-regressions.feature", import.meta.url);
  test("Control+click adds a branch to the filter, and Select works again after Filter (public#3298)", {tag: ["@viewers", "@realizes:charts.viewer.tree"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(16, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(17, "When user clicks on first \"Tree\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Tree\" card in \"Add Viewer\" dialog")));
    await session.step(18, "Then tree viewer should be visible", () => shouldBe(page, el("tree viewer"), "visible"));
    await session.step(19, "When user clicks on grid", () => clickOn(page, el("grid")));
    await session.step(20, "And user clicks on settings icon of tree viewer", () => clickOn(page, el("settings icon of tree viewer")));
    await session.step(21, "And user clicks on \"...\" button in \"Hierarchy\" property in context panel", () => clickOn(page, el("\"...\" button in \"Hierarchy\" property in context panel")));
    await session.step(22, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(23, "When user clicks on \"None\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"None\" link in \"Select columns...\" dialog")));
    await session.step(24, "Then \"0 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"0 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(25, "When user toggles the \"CONTROL\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "CONTROL", el("\"Select columns...\" dialog")));
    await session.step(26, "And user toggles the \"SEX\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "SEX", el("\"Select columns...\" dialog")));
    await session.step(27, "And user toggles the \"RACE\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "RACE", el("\"Select columns...\" dialog")));
    await session.step(28, "Then \"3 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"3 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(29, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(30, "Then \"Select columns...\" dialog should be absent", () => shouldBe(page, el("\"Select columns...\" dialog"), "absent"));
    await session.step(31, "And \"Hierarchy\" property of tree viewer should be \"CONTROL, SEX, RACE\"", () => propertyShouldBe(page, "Hierarchy", el("tree viewer"), "CONTROL, SEX, RACE"));
    await session.step(32, "And tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(35, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
    await session.step(36, "When user selects \"Filter\" in \"On Click\" property in context panel", () => selectIn(page, "Filter", el("\"On Click\" property in context panel")));
    await session.step(37, "Then \"On Click\" property of tree viewer should be \"Filter\"", () => propertyShouldBe(page, "On Click", el("tree viewer"), "Filter"));
    await session.step(38, "When user clicks on the \"branch All | false | F | Asian\" area of tree viewer", () => clickArea(page, "branch All | false | F | Asian", el("tree viewer")));
    await session.step(39, "Then 7 rows should pass the filter", () => filterPasses(page, 7));
    await session.step(40, "When user clicks on the \"branch All | false | M | Asian\" area of tree viewer holding Control", () => clickAreaHolding(page, "branch All | false | M | Asian", el("tree viewer"), "Control"));
    await session.step(41, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
    await session.step(42, "And the filter should pass exactly the rows where \"RACE\" is \"Asian\"", () => filterIsExactlyCategory(page, "RACE", "Asian"));
    await session.step(43, "And no rows should be selected", () => noneSelected(page));
    await session.step(44, "When user selects \"Select\" in \"On Click\" property in context panel", () => selectIn(page, "Select", el("\"On Click\" property in context panel")));
    await session.step(45, "Then \"On Click\" property of tree viewer should be \"Select\"", () => propertyShouldBe(page, "On Click", el("tree viewer"), "Select"));
    await session.step(46, "When user resets the filter", () => resetFilter(page));
    await session.step(47, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
    await session.step(48, "When user clicks on the \"branch All | false | F | Asian\" area of tree viewer", () => clickArea(page, "branch All | false | F | Asian", el("tree viewer")));
    await session.step(49, "Then 7 rows should be selected", () => selectedRowCount(page, 7));
    await session.step(50, "And no rows where \"SEX\" is \"M\" should be selected", () => noneOfSelected(page, "SEX", "M"));
    await session.step(51, "When user clicks on the \"branch All | false | M | Asian\" area of tree viewer holding Control", () => clickAreaHolding(page, "branch All | false | M | Asian", el("tree viewer"), "Control"));
    await session.step(52, "Then 15 rows should be selected", () => selectedRowCount(page, 15));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
  });
  test("Rows selected in the grid are drawn on the tree (public#3249)", {tag: ["@viewers", "@realizes:charts.viewer.tree"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(16, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(17, "When user clicks on first \"Tree\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Tree\" card in \"Add Viewer\" dialog")));
    await session.step(18, "Then tree viewer should be visible", () => shouldBe(page, el("tree viewer"), "visible"));
    await session.step(19, "When user clicks on grid", () => clickOn(page, el("grid")));
    await session.step(20, "And user clicks on settings icon of tree viewer", () => clickOn(page, el("settings icon of tree viewer")));
    await session.step(21, "And user clicks on \"...\" button in \"Hierarchy\" property in context panel", () => clickOn(page, el("\"...\" button in \"Hierarchy\" property in context panel")));
    await session.step(22, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(23, "When user clicks on \"None\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"None\" link in \"Select columns...\" dialog")));
    await session.step(24, "Then \"0 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"0 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(25, "When user toggles the \"CONTROL\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "CONTROL", el("\"Select columns...\" dialog")));
    await session.step(26, "And user toggles the \"SEX\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "SEX", el("\"Select columns...\" dialog")));
    await session.step(27, "And user toggles the \"RACE\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "RACE", el("\"Select columns...\" dialog")));
    await session.step(28, "Then \"3 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"3 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(29, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(30, "Then \"Select columns...\" dialog should be absent", () => shouldBe(page, el("\"Select columns...\" dialog"), "absent"));
    await session.step(31, "And \"Hierarchy\" property of tree viewer should be \"CONTROL, SEX, RACE\"", () => propertyShouldBe(page, "Hierarchy", el("tree viewer"), "CONTROL, SEX, RACE"));
    await session.step(32, "And tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(56, "When user takes a snapshot of tree viewer", () => takeSnapshot(page, el("tree viewer")));
    await session.step(57, "And user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
    await session.step(58, "Then 15 rows should be selected", () => selectedRowCount(page, 15));
    await session.step(59, "And tree viewer should show more selection highlight than before", () => moreHighlight(page, el("tree viewer")));
    await session.step(60, "When user takes a snapshot of tree viewer", () => takeSnapshot(page, el("tree viewer")));
    await session.step(61, "And user clears the row selection", () => clearSelection(page));
    await session.step(62, "Then no rows should be selected", () => noneSelected(page));
    await session.step(63, "And tree viewer should show less selection highlight than before", () => lessHighlight(page, el("tree viewer")));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
  });
});
