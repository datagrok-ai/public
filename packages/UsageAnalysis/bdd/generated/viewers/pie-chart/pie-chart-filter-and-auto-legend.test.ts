/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pie-chart/pie-chart-filter-and-auto-legend.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pie-chart, core#982]
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
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterIsExactlyCategory, filterPassesAll, noneOfFiltered} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, legendLists, noErrors, readingIs, setProperties, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pie chart's own filter and the Auto legend", () => {
  const session = feature(test, "features/viewers/pie-chart/pie-chart-filter-and-auto-legend.feature", import.meta.url);
  test("Pie chart's own filter and the Auto legend", {tag: ["@journey", "@viewers", "@realizes:viewers.pie-chart", "@realizes:core#982"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"]]), [["Category","RACE"]]);
    await session.step(16, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(17, "And pie chart viewer should show 1000 rows", () => showsRows(page, el("pie chart viewer"), 1000));
    await run.scenario("The pie's own Filter narrows the wedges and leaves the table filter alone, and a wedge click under it filters the table to that category alone (core#982)", async () => {
      await session.step(21, "Then the \"angle value of Caucasian\" reading of pie chart viewer should be 896", () => readingIs(page, "angle value of Caucasian", el("pie chart viewer"), 896));
      await session.step(22, "When user sets \"Filter\" property of pie chart viewer to \"${AGE} > 40\"", () => setProperty(page, "Filter", el("pie chart viewer"), "${AGE} > 40"));
      await session.step(23, "Then pie chart viewer should show 635 rows", () => showsRows(page, el("pie chart viewer"), 635));
      await session.step(24, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
      await session.step(25, "And the \"angle value of Caucasian\" reading of pie chart viewer should be 562", () => readingIs(page, "angle value of Caucasian", el("pie chart viewer"), 562));
      await session.step(26, "And the \"angle value of Other\" reading of pie chart viewer should be 48", () => readingIs(page, "angle value of Other", el("pie chart viewer"), 48));
      await session.step(27, "And the \"angle value of Black\" reading of pie chart viewer should be 19", () => readingIs(page, "angle value of Black", el("pie chart viewer"), 19));
      await session.step(28, "And the \"angle value of Asian\" reading of pie chart viewer should be 6", () => readingIs(page, "angle value of Asian", el("pie chart viewer"), 6));
      await session.step(29, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(30, "When user sets \"On Click\" property of pie chart viewer to \"Filter\"", () => setProperty(page, "On Click", el("pie chart viewer"), "Filter"));
      await session.step(31, "And user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
      await session.step(32, "Then the filter should pass exactly the rows where \"RACE\" is \"Caucasian\"", () => filterIsExactlyCategory(page, "RACE", "Caucasian"));
      await session.step(33, "And no rows where \"RACE\" is \"Asian\" should pass the filter", () => noneOfFiltered(page, "RACE", "Asian"));
      await session.step(34, "When user clicks on the \"empty space\" area of pie chart viewer", () => clickArea(page, "empty space", el("pie chart viewer")));
      await session.step(35, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(36, "When user sets \"On Click\" property of pie chart viewer to \"Select\"", () => setProperty(page, "On Click", el("pie chart viewer"), "Select"));
      await session.step(37, "And user sets \"Filter\" property of pie chart viewer to \"\"", () => setProperty(page, "Filter", el("pie chart viewer"), ""));
      await session.step(38, "Then pie chart viewer should show 1000 rows", () => showsRows(page, el("pie chart viewer"), 1000));
      await session.step(39, "And the \"angle value of Caucasian\" reading of pie chart viewer should be 896", () => readingIs(page, "angle value of Caucasian", el("pie chart viewer"), 896));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Legend Visibility Auto shows the legend while the wedges carry only percentages and hides it once they are labelled (GROK-20793)", async () => {
      await session.step(43, "When user sets properties of pie chart viewer:", () => setProperties(page, el("pie chart viewer"), [["Legend Visibility","Auto"],["Show Label","false"]]), [["Legend Visibility","Auto"],["Show Label","false"]]);
      await session.step(46, "Then legend of pie chart viewer should be visible", () => shouldBe(page, el("legend of pie chart viewer"), "visible"));
      await session.step(47, "And the legend of pie chart viewer should list 4 items", () => legendLists(page, el("pie chart viewer"), 4));
      await session.step(48, "When user sets \"Show Label\" property of pie chart viewer to \"true\"", () => setProperty(page, "Show Label", el("pie chart viewer"), "true"));
      await session.step(49, "Then legend of pie chart viewer should be hidden", () => shouldBe(page, el("legend of pie chart viewer"), "hidden"));
      await session.step(50, "When user sets \"Show Label\" property of pie chart viewer to \"false\"", () => setProperty(page, "Show Label", el("pie chart viewer"), "false"));
      await session.step(51, "Then legend of pie chart viewer should be visible", () => shouldBe(page, el("legend of pie chart viewer"), "visible"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
