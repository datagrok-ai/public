/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tree/tree-on-click-filter.feature
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
import {filterIsExactlyCategory, filterPasses, noneOfFiltered, noneSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, clickAreaHolding, hasArea, noErrors, propertyShouldBe} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Tree with On Click = Filter filters the table to the clicked branch", () => {
  const session = feature(test, "features/viewers/tree/tree-on-click-filter.feature", import.meta.url);
  test("A click filters to the branch, Shift+click adds a branch, a plain click starts over", {tag: ["@viewers", "@realizes:charts.viewer.tree"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(19, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(20, "When user clicks on first \"Tree\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Tree\" card in \"Add Viewer\" dialog")));
    await session.step(21, "Then tree viewer should be visible", () => shouldBe(page, el("tree viewer"), "visible"));
    await session.step(22, "When user clicks on grid", () => clickOn(page, el("grid")));
    await session.step(23, "And user clicks on settings icon of tree viewer", () => clickOn(page, el("settings icon of tree viewer")));
    await session.step(24, "And user clicks on \"...\" button in \"Hierarchy\" property in context panel", () => clickOn(page, el("\"...\" button in \"Hierarchy\" property in context panel")));
    await session.step(25, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(26, "When user clicks on \"None\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"None\" link in \"Select columns...\" dialog")));
    await session.step(27, "Then \"0 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"0 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(28, "When user toggles the \"CONTROL\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "CONTROL", el("\"Select columns...\" dialog")));
    await session.step(29, "And user toggles the \"SEX\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "SEX", el("\"Select columns...\" dialog")));
    await session.step(30, "And user toggles the \"RACE\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "RACE", el("\"Select columns...\" dialog")));
    await session.step(31, "Then \"3 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"3 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(32, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(33, "Then \"Select columns...\" dialog should be absent", () => shouldBe(page, el("\"Select columns...\" dialog"), "absent"));
    await session.step(34, "And \"Hierarchy\" property of tree viewer should be \"CONTROL, SEX, RACE\"", () => propertyShouldBe(page, "Hierarchy", el("tree viewer"), "CONTROL, SEX, RACE"));
    await session.step(35, "And tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(38, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
    await session.step(39, "When user selects \"Filter\" in \"On Click\" property in context panel", () => selectIn(page, "Filter", el("\"On Click\" property in context panel")));
    await session.step(40, "Then \"On Click\" property of tree viewer should be \"Filter\"", () => propertyShouldBe(page, "On Click", el("tree viewer"), "Filter"));
    await session.step(41, "And 1000 rows should pass the filter", () => filterPasses(page, 1000));
    await session.step(42, "When user clicks on the \"branch All | false | F | Asian\" area of tree viewer", () => clickArea(page, "branch All | false | F | Asian", el("tree viewer")));
    await session.step(43, "Then 7 rows should pass the filter", () => filterPasses(page, 7));
    await session.step(44, "And no rows where \"SEX\" is \"M\" should pass the filter", () => noneOfFiltered(page, "SEX", "M"));
    await session.step(45, "And no rows where \"RACE\" is \"Caucasian\" should pass the filter", () => noneOfFiltered(page, "RACE", "Caucasian"));
    await session.step(46, "And no rows should be selected", () => noneSelected(page));
    await session.step(47, "When user clicks on the \"branch All | false | M | Asian\" area of tree viewer holding Shift", () => clickAreaHolding(page, "branch All | false | M | Asian", el("tree viewer"), "Shift"));
    await session.step(48, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
    await session.step(49, "And the filter should pass exactly the rows where \"RACE\" is \"Asian\"", () => filterIsExactlyCategory(page, "RACE", "Asian"));
    await session.step(50, "When user clicks on the \"branch All | true | F | Caucasian\" area of tree viewer", () => clickArea(page, "branch All | true | F | Caucasian", el("tree viewer")));
    await session.step(51, "Then 5 rows should pass the filter", () => filterPasses(page, 5));
    await session.step(52, "And no rows where \"CONTROL\" is \"false\" should pass the filter", () => noneOfFiltered(page, "CONTROL", "false"));
    await session.step(53, "And no rows where \"RACE\" is \"Asian\" should pass the filter", () => noneOfFiltered(page, "RACE", "Asian"));
    await session.step(54, "And no rows should be selected", () => noneSelected(page));
    await session.step(55, "And no errors should have been logged", () => noErrors(page));
  });
  test("A first click made with Shift filters to the branch as a plain click does", {tag: ["@viewers", "@realizes:charts.viewer.tree", "@known-failure:66"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(19, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(20, "When user clicks on first \"Tree\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Tree\" card in \"Add Viewer\" dialog")));
    await session.step(21, "Then tree viewer should be visible", () => shouldBe(page, el("tree viewer"), "visible"));
    await session.step(22, "When user clicks on grid", () => clickOn(page, el("grid")));
    await session.step(23, "And user clicks on settings icon of tree viewer", () => clickOn(page, el("settings icon of tree viewer")));
    await session.step(24, "And user clicks on \"...\" button in \"Hierarchy\" property in context panel", () => clickOn(page, el("\"...\" button in \"Hierarchy\" property in context panel")));
    await session.step(25, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(26, "When user clicks on \"None\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"None\" link in \"Select columns...\" dialog")));
    await session.step(27, "Then \"0 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"0 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(28, "When user toggles the \"CONTROL\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "CONTROL", el("\"Select columns...\" dialog")));
    await session.step(29, "And user toggles the \"SEX\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "SEX", el("\"Select columns...\" dialog")));
    await session.step(30, "And user toggles the \"RACE\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "RACE", el("\"Select columns...\" dialog")));
    await session.step(31, "Then \"3 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"3 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(32, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(33, "Then \"Select columns...\" dialog should be absent", () => shouldBe(page, el("\"Select columns...\" dialog"), "absent"));
    await session.step(34, "And \"Hierarchy\" property of tree viewer should be \"CONTROL, SEX, RACE\"", () => propertyShouldBe(page, "Hierarchy", el("tree viewer"), "CONTROL, SEX, RACE"));
    await session.step(35, "And tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await knownFailure(async () => {
      await session.step(61, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
      await session.step(62, "When user selects \"Filter\" in \"On Click\" property in context panel", () => selectIn(page, "Filter", el("\"On Click\" property in context panel")));
      await session.step(63, "Then \"On Click\" property of tree viewer should be \"Filter\"", () => propertyShouldBe(page, "On Click", el("tree viewer"), "Filter"));
      await session.step(64, "And 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(65, "When user clicks on the \"branch All | false | F | Asian\" area of tree viewer holding Shift", () => clickAreaHolding(page, "branch All | false | F | Asian", el("tree viewer"), "Shift"));
      await session.step(66, "Then 7 rows should pass the filter", () => filterPasses(page, 7));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    }, 66);
  });
});
