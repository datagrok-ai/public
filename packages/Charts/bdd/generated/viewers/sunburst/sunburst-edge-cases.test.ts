/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/sunburst/sunburst-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.sunburst, public#4425]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, filterPassesAll, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, noErrors, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Sunburst edge cases: an empty hierarchy, a hierarchy column removed, On Click switched from Filter to Select", () => {
  const session = feature(test, "features/viewers/sunburst/sunburst-edge-cases.feature", import.meta.url);
  test("A hierarchy with no column shows the viewer's message instead of a chart", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.sunburst"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a sunburst viewer with:", () => addViewerWith(page, "sunburst", [["hierarchyColumnNames","SEX, RACE"]]), [["hierarchyColumnNames","SEX, RACE"]]);
    await session.step(19, "Then the \"hierarchy columns\" reading of sunburst viewer should be \"SEX, RACE\"", () => readingReads(page, "hierarchy columns", el("sunburst viewer"), "SEX, RACE"));
    await session.step(20, "And the \"segments\" reading of sunburst viewer should be 10", () => readingIs(page, "segments", el("sunburst viewer"), 10));
    await session.step(23, "When user sets \"hierarchyColumnNames\" property of sunburst viewer to \"\"", () => setProperty(page, "hierarchyColumnNames", el("sunburst viewer"), ""));
    await session.step(24, "Then \"The Sunburst viewer requires at least one categorical column with fewer than 500 unique categories\" text in sunburst viewer should be visible", () => shouldBe(page, el("\"The Sunburst viewer requires at least one categorical column with fewer than 500 unique categories\" text in sunburst viewer"), "visible"));
    await session.step(25, "When user sets \"hierarchyColumnNames\" property of sunburst viewer to \"SEX, RACE\"", () => setProperty(page, "hierarchyColumnNames", el("sunburst viewer"), "SEX, RACE"));
    await session.step(26, "Then \"The Sunburst viewer requires at least one categorical column with fewer than 500 unique categories\" text in sunburst viewer should be absent", () => shouldBe(page, el("\"The Sunburst viewer requires at least one categorical column with fewer than 500 unique categories\" text in sunburst viewer"), "absent"));
    await session.step(27, "And the \"segments\" reading of sunburst viewer should be 10", () => readingIs(page, "segments", el("sunburst viewer"), 10));
    await session.step(28, "And no errors should have been logged", () => noErrors(page));
  });
  test("A hierarchy column removed from the table is dropped from the hierarchy", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.sunburst"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a sunburst viewer with:", () => addViewerWith(page, "sunburst", [["hierarchyColumnNames","SEX, RACE"]]), [["hierarchyColumnNames","SEX, RACE"]]);
    await session.step(19, "Then the \"hierarchy columns\" reading of sunburst viewer should be \"SEX, RACE\"", () => readingReads(page, "hierarchy columns", el("sunburst viewer"), "SEX, RACE"));
    await session.step(20, "And the \"segments\" reading of sunburst viewer should be 10", () => readingIs(page, "segments", el("sunburst viewer"), 10));
    await session.step(31, "When user removes \"SEX\" column", () => removeColumn(page, "SEX"));
    await session.step(32, "Then the \"hierarchy columns\" reading of sunburst viewer should be \"RACE\"", () => readingReads(page, "hierarchy columns", el("sunburst viewer"), "RACE"));
    await session.step(33, "And the \"segments\" reading of sunburst viewer should be 4", () => readingIs(page, "segments", el("sunburst viewer"), 4));
    await session.step(34, "And the \"rows of segment Caucasian\" reading of sunburst viewer should be 896", () => readingIs(page, "rows of segment Caucasian", el("sunburst viewer"), 896));
    await session.step(35, "And no errors should have been logged", () => noErrors(page));
  });
  test("Switching On Click from Filter to Select gives back the rows its clicked segment filtered out", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.sunburst", "@known-failure:43", "@realizes:public#4425"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a sunburst viewer with:", () => addViewerWith(page, "sunburst", [["hierarchyColumnNames","SEX, RACE"]]), [["hierarchyColumnNames","SEX, RACE"]]);
    await session.step(19, "Then the \"hierarchy columns\" reading of sunburst viewer should be \"SEX, RACE\"", () => readingReads(page, "hierarchy columns", el("sunburst viewer"), "SEX, RACE"));
    await session.step(20, "And the \"segments\" reading of sunburst viewer should be 10", () => readingIs(page, "segments", el("sunburst viewer"), 10));
    await knownFailure(async () => {
      await session.step(39, "When user sets \"onClick\" property of sunburst viewer to \"Filter\"", () => setProperty(page, "onClick", el("sunburst viewer"), "Filter"));
      await session.step(40, "And user clicks on the \"segment F\" area of sunburst viewer", () => clickArea(page, "segment F", el("sunburst viewer")));
      await session.step(41, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(42, "When user sets \"onClick\" property of sunburst viewer to \"Select\"", () => setProperty(page, "onClick", el("sunburst viewer"), "Select"));
      await session.step(43, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    }, 43);
  });
});
