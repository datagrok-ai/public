/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tree/tree-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.tree, public#4425]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {filterPasses, filterPassesAll, filterTo, removeColumn, renameColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, hasArea, hasNoArea, noErrors, propertyShouldBe, readingIs, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Tree edge cases: a column renamed under a filter, a hierarchy column removed, empty values, On Click switched from Filter to Select", () => {
  const session = feature(test, "features/viewers/tree/tree-edge-cases.feature", import.meta.url);
  test("A hierarchy column renamed under a filter of fewer than 500 rows is followed", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.tree"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","CONTROL, SEX, RACE"]]), [["hierarchyColumnNames","CONTROL, SEX, RACE"]]);
    await session.step(20, "Then tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(21, "And the \"rows shown\" reading of tree viewer should be 1000", () => readingIs(page, "rows shown", el("tree viewer"), 1000));
    await session.step(24, "When user filters rows where \"SEX\" is \"M\"", () => filterTo(page, "SEX", "M"));
    await session.step(25, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
    await session.step(26, "When user renames \"RACE\" column to \"ORIGIN\"", () => renameColumn(page, "RACE", "ORIGIN"));
    await session.step(27, "Then \"hierarchyColumnNames\" property of tree viewer should be \"CONTROL, SEX, ORIGIN\"", () => propertyShouldBe(page, "hierarchyColumnNames", el("tree viewer"), "CONTROL, SEX, ORIGIN"));
    await session.step(28, "And tree viewer should have a \"branch All | false | M | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | M | Asian"));
    await session.step(29, "And the \"rows shown\" reading of tree viewer should be 447", () => readingIs(page, "rows shown", el("tree viewer"), 447));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
  });
  test("A hierarchy column removed from the table is dropped from the hierarchy", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.tree"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","CONTROL, SEX, RACE"]]), [["hierarchyColumnNames","CONTROL, SEX, RACE"]]);
    await session.step(20, "Then tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(21, "And the \"rows shown\" reading of tree viewer should be 1000", () => readingIs(page, "rows shown", el("tree viewer"), 1000));
    await session.step(33, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
    await session.step(34, "Then \"hierarchyColumnNames\" property of tree viewer should be \"CONTROL, SEX\"", () => propertyShouldBe(page, "hierarchyColumnNames", el("tree viewer"), "CONTROL, SEX"));
    await session.step(35, "And tree viewer should have a \"branch All | false | F\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F"));
    await session.step(36, "And tree viewer should not have a \"branch All | false | F | Asian\" area", () => hasNoArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
  test("An empty value is the (empty) branch only while Include Nulls is on", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.tree"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","CONTROL, SEX, RACE"]]), [["hierarchyColumnNames","CONTROL, SEX, RACE"]]);
    await session.step(20, "Then tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(21, "And the \"rows shown\" reading of tree viewer should be 1000", () => readingIs(page, "rows shown", el("tree viewer"), 1000));
    await session.step(40, "Given user opens a table \"with-empty\" with:", () => openTableOf(page, "with-empty", [["site","arm"],["A","x"],["A",""],["B","y"]]), [["site","arm"],["A","x"],["A",""],["B","y"]]);
    await session.step(45, "And user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","site, arm"],["includeNulls","true"]]), [["hierarchyColumnNames","site, arm"],["includeNulls","true"]]);
    await session.step(48, "Then tree viewer should have a \"branch All | A | (empty)\" area", () => hasArea(page, el("tree viewer"), "branch All | A | (empty)"));
    await session.step(49, "When user sets \"includeNulls\" property of tree viewer to \"false\"", () => setProperty(page, "includeNulls", el("tree viewer"), "false"));
    await session.step(50, "Then tree viewer should not have a \"branch All | A | (empty)\" area", () => hasNoArea(page, el("tree viewer"), "branch All | A | (empty)"));
    await session.step(51, "And tree viewer should have a \"branch All | A | x\" area", () => hasArea(page, el("tree viewer"), "branch All | A | x"));
    await session.step(52, "When user sets \"includeNulls\" property of tree viewer to \"true\"", () => setProperty(page, "includeNulls", el("tree viewer"), "true"));
    await session.step(53, "Then tree viewer should have a \"branch All | A | (empty)\" area", () => hasArea(page, el("tree viewer"), "branch All | A | (empty)"));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
  });
  test("Switching On Click from Filter to Select gives back the rows its clicked branch filtered out", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.tree", "@known-failure:62", "@realizes:public#4425"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","CONTROL, SEX, RACE"]]), [["hierarchyColumnNames","CONTROL, SEX, RACE"]]);
    await session.step(20, "Then tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(21, "And the \"rows shown\" reading of tree viewer should be 1000", () => readingIs(page, "rows shown", el("tree viewer"), 1000));
    await knownFailure(async () => {
      await session.step(58, "When user sets \"onClick\" property of tree viewer to \"Filter\"", () => setProperty(page, "onClick", el("tree viewer"), "Filter"));
      await session.step(59, "And user clicks on the \"branch All | false | F | Asian\" area of tree viewer", () => clickArea(page, "branch All | false | F | Asian", el("tree viewer")));
      await session.step(60, "Then 7 rows should pass the filter", () => filterPasses(page, 7));
      await session.step(61, "When user sets \"onClick\" property of tree viewer to \"Select\"", () => setProperty(page, "onClick", el("tree viewer"), "Select"));
      await session.step(62, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    }, 62);
  });
});
