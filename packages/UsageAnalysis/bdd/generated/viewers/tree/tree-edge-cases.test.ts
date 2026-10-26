/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tree/tree-edge-cases.feature
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
import {clickOn, enterInto, finishedUpdating, isExpanded} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterBetween, filterPasses, filterTo, removeColumn, renameColumn, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, hasNoArea, noErrors, propertyShouldBe, readingIs, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Tree edge cases: renamed and removed columns, empty values, no rows, out-of-range style values", () => {
  const session = feature(test, "features/viewers/tree/tree-edge-cases.feature", import.meta.url);
  test("A hierarchy column renamed under a filter of fewer than 500 rows is followed", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.tree"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","CONTROL, SEX, RACE"]]), [["hierarchyColumnNames","CONTROL, SEX, RACE"]]);
    await session.step(22, "Then tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(23, "And the \"rows shown\" reading of tree viewer should be 1000", () => readingIs(page, "rows shown", el("tree viewer"), 1000));
    await session.step(26, "When user filters rows where \"SEX\" is \"M\"", () => filterTo(page, "SEX", "M"));
    await session.step(27, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
    await session.step(28, "When user renames \"RACE\" column to \"ORIGIN\"", () => renameColumn(page, "RACE", "ORIGIN"));
    await session.step(29, "Then \"hierarchyColumnNames\" property of tree viewer should be \"CONTROL, SEX, ORIGIN\"", () => propertyShouldBe(page, "hierarchyColumnNames", el("tree viewer"), "CONTROL, SEX, ORIGIN"));
    await session.step(30, "And tree viewer should have a \"branch All | false | M | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | M | Asian"));
    await session.step(31, "And the \"rows shown\" reading of tree viewer should be 447", () => readingIs(page, "rows shown", el("tree viewer"), 447));
    await session.step(32, "And no errors should have been logged", () => noErrors(page));
  });
  test("A hierarchy column renamed with every row passing is followed", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.tree"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","CONTROL, SEX, RACE"]]), [["hierarchyColumnNames","CONTROL, SEX, RACE"]]);
    await session.step(22, "Then tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(23, "And the \"rows shown\" reading of tree viewer should be 1000", () => readingIs(page, "rows shown", el("tree viewer"), 1000));
    await session.step(35, "When user renames \"RACE\" column to \"ORIGIN\"", () => renameColumn(page, "RACE", "ORIGIN"));
    await session.step(36, "Then \"hierarchyColumnNames\" property of tree viewer should be \"CONTROL, SEX, ORIGIN\"", () => propertyShouldBe(page, "hierarchyColumnNames", el("tree viewer"), "CONTROL, SEX, ORIGIN"));
    await session.step(37, "And tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
  test("A hierarchy column removed from the table is dropped from the hierarchy", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.tree"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","CONTROL, SEX, RACE"]]), [["hierarchyColumnNames","CONTROL, SEX, RACE"]]);
    await session.step(22, "Then tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(23, "And the \"rows shown\" reading of tree viewer should be 1000", () => readingIs(page, "rows shown", el("tree viewer"), 1000));
    await session.step(41, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
    await session.step(42, "Then \"hierarchyColumnNames\" property of tree viewer should be \"CONTROL, SEX\"", () => propertyShouldBe(page, "hierarchyColumnNames", el("tree viewer"), "CONTROL, SEX"));
    await session.step(43, "And tree viewer should have a \"branch All | false | F\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F"));
    await session.step(44, "And tree viewer should not have a \"branch All | false | F | Asian\" area", () => hasNoArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
  });
  test("Removing the Size column logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.tree", "@known-failure:55", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","CONTROL, SEX, RACE"]]), [["hierarchyColumnNames","CONTROL, SEX, RACE"]]);
    await session.step(22, "Then tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(23, "And the \"rows shown\" reading of tree viewer should be 1000", () => readingIs(page, "rows shown", el("tree viewer"), 1000));
    await knownFailure(async () => {
      await session.step(52, "When user sets \"sizeColumnName\" property of tree viewer to \"HEIGHT\"", () => setProperty(page, "sizeColumnName", el("tree viewer"), "HEIGHT"));
      await session.step(53, "And user removes \"HEIGHT\" column", () => removeColumn(page, "HEIGHT"));
      await session.step(54, "Then the table should not have a column \"HEIGHT\"", () => hasNoColumn(page, "HEIGHT"));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
      await session.step(56, "And tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    }, 55);
  });
  test("Removing the Color column logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.tree", "@known-failure:66", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","CONTROL, SEX, RACE"]]), [["hierarchyColumnNames","CONTROL, SEX, RACE"]]);
    await session.step(22, "Then tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(23, "And the \"rows shown\" reading of tree viewer should be 1000", () => readingIs(page, "rows shown", el("tree viewer"), 1000));
    await knownFailure(async () => {
      await session.step(63, "When user sets \"colorColumnName\" property of tree viewer to \"AGE\"", () => setProperty(page, "colorColumnName", el("tree viewer"), "AGE"));
      await session.step(64, "And user removes \"AGE\" column", () => removeColumn(page, "AGE"));
      await session.step(65, "Then the table should not have a column \"AGE\"", () => hasNoColumn(page, "AGE"));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
      await session.step(67, "And tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    }, 66);
  });
  test("A filter that removes every row logs no error and the branches come back", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.tree"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","CONTROL, SEX, RACE"]]), [["hierarchyColumnNames","CONTROL, SEX, RACE"]]);
    await session.step(22, "Then tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(23, "And the \"rows shown\" reading of tree viewer should be 1000", () => readingIs(page, "rows shown", el("tree viewer"), 1000));
    await session.step(70, "When user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(71, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(72, "And the \"rows shown\" reading of tree viewer should be 0", () => readingIs(page, "rows shown", el("tree viewer"), 0));
    await session.step(73, "And tree viewer should not have a \"branch All | false | F | Asian\" area", () => hasNoArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(74, "When user resets the filter", () => resetFilter(page));
    await session.step(75, "Then the \"rows shown\" reading of tree viewer should be 1000", () => readingIs(page, "rows shown", el("tree viewer"), 1000));
    await session.step(76, "And tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(77, "And no errors should have been logged", () => noErrors(page));
  });
  test("An empty value is the (empty) branch only while Include Nulls is on", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.tree"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","CONTROL, SEX, RACE"]]), [["hierarchyColumnNames","CONTROL, SEX, RACE"]]);
    await session.step(22, "Then tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(23, "And the \"rows shown\" reading of tree viewer should be 1000", () => readingIs(page, "rows shown", el("tree viewer"), 1000));
    await session.step(80, "Given user opens a table \"with-empty\" with:", () => openTableOf(page, "with-empty", [["site","arm"],["A","x"],["A",""],["B","y"]]), [["site","arm"],["A","x"],["A",""],["B","y"]]);
    await session.step(85, "And user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","site, arm"],["includeNulls","true"]]), [["hierarchyColumnNames","site, arm"],["includeNulls","true"]]);
    await session.step(88, "Then tree viewer should have a \"branch All | A | (empty)\" area", () => hasArea(page, el("tree viewer"), "branch All | A | (empty)"));
    await session.step(89, "When user sets \"includeNulls\" property of tree viewer to \"false\"", () => setProperty(page, "includeNulls", el("tree viewer"), "false"));
    await session.step(90, "Then tree viewer should not have a \"branch All | A | (empty)\" area", () => hasNoArea(page, el("tree viewer"), "branch All | A | (empty)"));
    await session.step(91, "And tree viewer should have a \"branch All | A | x\" area", () => hasArea(page, el("tree viewer"), "branch All | A | x"));
    await session.step(92, "When user sets \"includeNulls\" property of tree viewer to \"true\"", () => setProperty(page, "includeNulls", el("tree viewer"), "true"));
    await session.step(93, "Then tree viewer should have a \"branch All | A | (empty)\" area", () => hasArea(page, el("tree viewer"), "branch All | A | (empty)"));
    await session.step(94, "And no errors should have been logged", () => noErrors(page));
  });
  test("Values above the declared maximum typed in the context panel log no error", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.tree"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","CONTROL, SEX, RACE"]]), [["hierarchyColumnNames","CONTROL, SEX, RACE"]]);
    await session.step(22, "Then tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(23, "And the \"rows shown\" reading of tree viewer should be 1000", () => readingIs(page, "rows shown", el("tree viewer"), 1000));
    await session.step(97, "When user clicks on grid", () => clickOn(page, el("grid")));
    await session.step(98, "And user clicks on settings icon of tree viewer", () => clickOn(page, el("settings icon of tree viewer")));
    await session.step(99, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
    await session.step(100, "When user enters \"6\" into \"Initial Tree Depth\" property in context panel", () => enterInto(page, "6", el("\"Initial Tree Depth\" property in context panel")));
    await session.step(101, "And user enters \"31\" into \"Font Size\" property in context panel", () => enterInto(page, "31", el("\"Font Size\" property in context panel")));
    await session.step(102, "Then tree viewer should have finished updating", () => finishedUpdating(page, el("tree viewer")));
    await session.step(103, "And tree viewer should have a \"branch All | false | F | Asian\" area", () => hasArea(page, el("tree viewer"), "branch All | false | F | Asian"));
    await session.step(104, "And no errors should have been logged", () => noErrors(page));
  });
});
