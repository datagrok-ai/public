/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/sunburst/sunburst-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.sunburst]
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
import {filterBetween, filterPasses, filterPassesAll, removeColumn, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, hasNoArea, loadLayout, noErrors, readingIs, readingReads, saveLayout, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {doubleClickEmptySpace} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Sunburst edge cases: no rows, an empty hierarchy, a hierarchy column removed, a layout with a missing column", () => {
  const session = feature(test, "features/viewers/sunburst/sunburst-edge-cases.feature", import.meta.url);
  test("A filter that removes every row leaves the sunburst empty without an error", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.sunburst"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a sunburst viewer with:", () => addViewerWith(page, "sunburst", [["hierarchyColumnNames","SEX, RACE"]]), [["hierarchyColumnNames","SEX, RACE"]]);
    await session.step(22, "Then the \"hierarchy columns\" reading of sunburst viewer should be \"SEX, RACE\"", () => readingReads(page, "hierarchy columns", el("sunburst viewer"), "SEX, RACE"));
    await session.step(23, "And the \"segments\" reading of sunburst viewer should be 10", () => readingIs(page, "segments", el("sunburst viewer"), 10));
    await session.step(26, "When user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(27, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(28, "And the \"rows shown\" reading of sunburst viewer should be 0", () => readingIs(page, "rows shown", el("sunburst viewer"), 0));
    await session.step(29, "And sunburst viewer should not have a \"segment F\" area", () => hasNoArea(page, el("sunburst viewer"), "segment F"));
    await session.step(30, "When user resets the filter", () => resetFilter(page));
    await session.step(31, "Then the \"rows shown\" reading of sunburst viewer should be 1000", () => readingIs(page, "rows shown", el("sunburst viewer"), 1000));
    await session.step(32, "And the \"segments\" reading of sunburst viewer should be 10", () => readingIs(page, "segments", el("sunburst viewer"), 10));
    await session.step(33, "And no errors should have been logged", () => noErrors(page));
  });
  test("A hierarchy with no column shows the viewer's message instead of a chart", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.sunburst"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a sunburst viewer with:", () => addViewerWith(page, "sunburst", [["hierarchyColumnNames","SEX, RACE"]]), [["hierarchyColumnNames","SEX, RACE"]]);
    await session.step(22, "Then the \"hierarchy columns\" reading of sunburst viewer should be \"SEX, RACE\"", () => readingReads(page, "hierarchy columns", el("sunburst viewer"), "SEX, RACE"));
    await session.step(23, "And the \"segments\" reading of sunburst viewer should be 10", () => readingIs(page, "segments", el("sunburst viewer"), 10));
    await session.step(36, "When user sets \"hierarchyColumnNames\" property of sunburst viewer to \"\"", () => setProperty(page, "hierarchyColumnNames", el("sunburst viewer"), ""));
    await session.step(37, "Then \"The Sunburst viewer requires at least one categorical column with fewer than 500 unique categories\" text in sunburst viewer should be visible", () => shouldBe(page, el("\"The Sunburst viewer requires at least one categorical column with fewer than 500 unique categories\" text in sunburst viewer"), "visible"));
    await session.step(38, "When user sets \"hierarchyColumnNames\" property of sunburst viewer to \"SEX, RACE\"", () => setProperty(page, "hierarchyColumnNames", el("sunburst viewer"), "SEX, RACE"));
    await session.step(39, "Then \"The Sunburst viewer requires at least one categorical column with fewer than 500 unique categories\" text in sunburst viewer should be absent", () => shouldBe(page, el("\"The Sunburst viewer requires at least one categorical column with fewer than 500 unique categories\" text in sunburst viewer"), "absent"));
    await session.step(40, "And the \"segments\" reading of sunburst viewer should be 10", () => readingIs(page, "segments", el("sunburst viewer"), 10));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
  });
  test("A hierarchy column removed from the table is dropped from the hierarchy", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.sunburst"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a sunburst viewer with:", () => addViewerWith(page, "sunburst", [["hierarchyColumnNames","SEX, RACE"]]), [["hierarchyColumnNames","SEX, RACE"]]);
    await session.step(22, "Then the \"hierarchy columns\" reading of sunburst viewer should be \"SEX, RACE\"", () => readingReads(page, "hierarchy columns", el("sunburst viewer"), "SEX, RACE"));
    await session.step(23, "And the \"segments\" reading of sunburst viewer should be 10", () => readingIs(page, "segments", el("sunburst viewer"), 10));
    await session.step(44, "When user removes \"SEX\" column", () => removeColumn(page, "SEX"));
    await session.step(45, "Then the \"hierarchy columns\" reading of sunburst viewer should be \"RACE\"", () => readingReads(page, "hierarchy columns", el("sunburst viewer"), "RACE"));
    await session.step(46, "And the \"segments\" reading of sunburst viewer should be 4", () => readingIs(page, "segments", el("sunburst viewer"), 4));
    await session.step(47, "And the \"rows of segment Caucasian\" reading of sunburst viewer should be 896", () => readingIs(page, "rows of segment Caucasian", el("sunburst viewer"), 896));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
  });
  test("A double-click on empty space clears the filter after the clicked segment's column is removed", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.sunburst"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a sunburst viewer with:", () => addViewerWith(page, "sunburst", [["hierarchyColumnNames","SEX, RACE"]]), [["hierarchyColumnNames","SEX, RACE"]]);
    await session.step(22, "Then the \"hierarchy columns\" reading of sunburst viewer should be \"SEX, RACE\"", () => readingReads(page, "hierarchy columns", el("sunburst viewer"), "SEX, RACE"));
    await session.step(23, "And the \"segments\" reading of sunburst viewer should be 10", () => readingIs(page, "segments", el("sunburst viewer"), 10));
    await session.step(51, "When user sets \"onClick\" property of sunburst viewer to \"Filter\"", () => setProperty(page, "onClick", el("sunburst viewer"), "Filter"));
    await session.step(52, "And user clicks on the \"segment F\" area of sunburst viewer", () => clickArea(page, "segment F", el("sunburst viewer")));
    await session.step(53, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
    await session.step(54, "When user removes \"SEX\" column", () => removeColumn(page, "SEX"));
    await session.step(55, "Then the \"hierarchy columns\" reading of sunburst viewer should be \"RACE\"", () => readingReads(page, "hierarchy columns", el("sunburst viewer"), "RACE"));
    await session.step(56, "When user double-clicks on empty plot space of sunburst viewer", () => doubleClickEmptySpace(page, el("sunburst viewer")));
    await session.step(57, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
    await session.step(58, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout naming a hierarchy column the table no longer has is applied without an error", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.sunburst"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a sunburst viewer with:", () => addViewerWith(page, "sunburst", [["hierarchyColumnNames","SEX, RACE"]]), [["hierarchyColumnNames","SEX, RACE"]]);
    await session.step(22, "Then the \"hierarchy columns\" reading of sunburst viewer should be \"SEX, RACE\"", () => readingReads(page, "hierarchy columns", el("sunburst viewer"), "SEX, RACE"));
    await session.step(23, "And the \"segments\" reading of sunburst viewer should be 10", () => readingIs(page, "segments", el("sunburst viewer"), 10));
    await session.step(61, "When user saves the layout of the current table view", () => saveLayout(page));
    await session.step(62, "And user removes \"SEX\" column", () => removeColumn(page, "SEX"));
    await session.step(63, "And user loads the saved layout", () => loadLayout(page));
    await session.step(64, "Then sunburst viewer should be visible", () => shouldBe(page, el("sunburst viewer"), "visible"));
    await session.step(65, "And the \"hierarchy columns\" reading of sunburst viewer should be \"RACE\"", () => readingReads(page, "hierarchy columns", el("sunburst viewer"), "RACE"));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
  });
  test("Closing the sunburst gives back the rows its clicked segment filtered out", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.sunburst", "@known-failure:77", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a sunburst viewer with:", () => addViewerWith(page, "sunburst", [["hierarchyColumnNames","SEX, RACE"]]), [["hierarchyColumnNames","SEX, RACE"]]);
    await session.step(22, "Then the \"hierarchy columns\" reading of sunburst viewer should be \"SEX, RACE\"", () => readingReads(page, "hierarchy columns", el("sunburst viewer"), "SEX, RACE"));
    await session.step(23, "And the \"segments\" reading of sunburst viewer should be 10", () => readingIs(page, "segments", el("sunburst viewer"), 10));
    await knownFailure(async () => {
      await session.step(72, "When user sets \"onClick\" property of sunburst viewer to \"Filter\"", () => setProperty(page, "onClick", el("sunburst viewer"), "Filter"));
      await session.step(73, "And user clicks on the \"segment F\" area of sunburst viewer", () => clickArea(page, "segment F", el("sunburst viewer")));
      await session.step(74, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(75, "When user clicks on close icon of sunburst viewer", () => clickOn(page, el("close icon of sunburst viewer")));
      await session.step(76, "Then the open tableview should have 0 sunburst viewers", () => viewerCount(page, 0, "sunburst"));
      await session.step(77, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
    }, 77);
  });
});
