/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/box-plot/box-plot-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot, core#307]
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
import {addCalculated, filterBetween, filterIsExactly, filterPassesAll, removeColumn, resetFilter, rowCount, tableColumnIncomplete} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, doubleClickArea, hasArea, hasNoArea, narrowerRange, noErrors, propertyShouldBe, rememberRange, rememberedRange, repainted, sameRange, scaleNarrower, scaleWider, setProperty, showsFewerRows, showsRows, widerRange} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaLies, zoomValueAxis} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Box plot filter semantics", () => {
  const session = feature(test, "features/viewers/box-plot/box-plot-filter.feature", import.meta.url);
  test("Box plot filter semantics", {tag: ["@journey", "@viewers", "@realizes:viewers.box-plot", "@realizes:core#307"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(12, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","Average Mass"],["Category 1","Series"]]), [["Value","Average Mass"],["Category 1","Series"]]);
    await session.step(15, "Then the table should have 100 rows", () => rowCount(page, 100));
    await run.scenario("The value range follows the filter", async () => {
      await session.step(18, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(19, "And \"Zoom Values By Filter\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Zoom Values By Filter", el("box plot viewer"), "true"));
      await session.step(20, "And box plot viewer should show 100 rows", () => showsRows(page, el("box plot viewer"), 100));
      await session.step(21, "When user filters rows where \"Average Mass\" is between 300 and 400", () => filterBetween(page, "Average Mass", 300, 400));
      await session.step(22, "Then the filter should pass exactly the rows where \"Average Mass\" is between 300 and 400", () => filterIsExactly(page, "Average Mass", 300, 400));
      await session.step(23, "And box plot viewer should show fewer rows than before", () => showsFewerRows(page, el("box plot viewer")));
      await session.step(24, "And box plot viewer should show a narrower value range than before", () => narrowerRange(page, el("box plot viewer")));
      await session.step(25, "When user resets the filter", () => resetFilter(page));
      await session.step(26, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(27, "And box plot viewer should show a wider value range than before", () => widerRange(page, el("box plot viewer")));
    });
    await run.scenario("Zoom Values By Filter off keeps the range", async () => {
      await session.step(30, "When user sets \"Zoom Values By Filter\" property of box plot viewer to \"false\"", () => setProperty(page, "Zoom Values By Filter", el("box plot viewer"), "false"));
      await session.step(31, "And user filters rows where \"Average Mass\" is between 300 and 400", () => filterBetween(page, "Average Mass", 300, 400));
      await session.step(32, "Then the filter should pass exactly the rows where \"Average Mass\" is between 300 and 400", () => filterIsExactly(page, "Average Mass", 300, 400));
      await session.step(33, "And box plot viewer should show the same value range as before", () => sameRange(page, el("box plot viewer")));
      await session.step(34, "When user sets \"Zoom Values By Filter\" property of box plot viewer to \"true\"", () => setProperty(page, "Zoom Values By Filter", el("box plot viewer"), "true"));
      await session.step(35, "Then box plot viewer should show a narrower value range than before", () => narrowerRange(page, el("box plot viewer")));
      await session.step(36, "When user resets the filter", () => resetFilter(page));
    });
    await run.scenario("The viewer's own filter leaves the table's alone", async () => {
      await session.step(39, "When user filters rows where \"Average Mass\" is between 300 and 400", () => filterBetween(page, "Average Mass", 300, 400));
      await session.step(40, "And user sets \"Filter\" property of box plot viewer to \"${Average Mass} > 350\"", () => setProperty(page, "Filter", el("box plot viewer"), "${Average Mass} > 350"));
      await session.step(41, "Then box plot viewer should show fewer rows than before", () => showsFewerRows(page, el("box plot viewer")));
      await session.step(42, "And the filter should pass exactly the rows where \"Average Mass\" is between 300 and 400", () => filterIsExactly(page, "Average Mass", 300, 400));
      await session.step(43, "When user sets \"Filter\" property of box plot viewer to \"\"", () => setProperty(page, "Filter", el("box plot viewer"), ""));
      await session.step(44, "And user resets the filter", () => resetFilter(page));
    });
    await run.scenario("Show Empty Categories drops and restores an empty-valued category", async () => {
      await session.step(47, "When user adds a calculated column \"AverageMassFixture\" with formula \"if(${Series} == \\\"Triazoles\\\", null, ${Average Mass})\"", () => addCalculated(page, "AverageMassFixture", "if(${Series} == \"Triazoles\", null, ${Average Mass})"));
      await session.step(48, "Then table \"spgi-100\" should have missing values in \"AverageMassFixture\" column", () => tableColumnIncomplete(page, "spgi-100", "AverageMassFixture"));
      await session.step(49, "When user sets \"Value\" property of box plot viewer to \"AverageMassFixture\"", () => setProperty(page, "Value", el("box plot viewer"), "AverageMassFixture"));
      await session.step(50, "And user sets \"Show Empty Categories\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Empty Categories", el("box plot viewer"), "true"));
      await session.step(51, "Then box plot viewer should have a \"category Triazoles\" area", () => hasArea(page, el("box plot viewer"), "category Triazoles"));
      await session.step(52, "When user sets \"Show Empty Categories\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Empty Categories", el("box plot viewer"), "false"));
      await session.step(53, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(54, "And box plot viewer should not have a \"category Triazoles\" area", () => hasNoArea(page, el("box plot viewer"), "category Triazoles"));
      await session.step(55, "When user sets \"Show Empty Categories\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Empty Categories", el("box plot viewer"), "true"));
      await session.step(56, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(57, "And box plot viewer should have a \"category Triazoles\" area", () => hasArea(page, el("box plot viewer"), "category Triazoles"));
      await session.step(58, "When user sets \"Value\" property of box plot viewer to \"Average Mass\"", () => setProperty(page, "Value", el("box plot viewer"), "Average Mass"));
      await session.step(59, "And user removes \"AverageMassFixture\" column", () => removeColumn(page, "AverageMassFixture"));
    });
    await run.scenario("A coloring under a filter keeps its color scale", async () => {
      await session.step(62, "When user filters rows where \"Average Mass\" is between 300 and 400", () => filterBetween(page, "Average Mass", 300, 400));
      await session.step(63, "And user sets \"Marker Color Column\" property of box plot viewer to \"TPSA\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "TPSA"));
      await session.step(64, "Then \"Marker Color Column\" property of box plot viewer should be \"TPSA\"", () => propertyShouldBe(page, "Marker Color Column", el("box plot viewer"), "TPSA"));
      await session.step(65, "And box plot viewer should have a \"color scale\" area", () => hasArea(page, el("box plot viewer"), "color scale"));
      await session.step(66, "And the \"color scale\" area of box plot viewer should lie outside the \"y axis\" area", () => areaLies(page, "color scale", el("box plot viewer"), "outside", "y axis"));
      await session.step(67, "And the \"color scale\" area of box plot viewer should lie outside the \"x axis\" area", () => areaLies(page, "color scale", el("box plot viewer"), "outside", "x axis"));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
      await session.step(69, "When user resets the filter", () => resetFilter(page));
      await session.step(70, "Then the color scale of box plot viewer should cover a wider range than before", () => scaleWider(page, el("box plot viewer")));
      await session.step(71, "When user filters rows where \"Average Mass\" is between 300 and 400", () => filterBetween(page, "Average Mass", 300, 400));
      await session.step(72, "Then the color scale of box plot viewer should cover a narrower range than before", () => scaleNarrower(page, el("box plot viewer")));
      await session.step(73, "When user resets the filter", () => resetFilter(page));
      await session.step(74, "And user sets \"Marker Color Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), ""));
    });
    await run.scenario("Toggling Show Empty Categories keeps the zoom (core#307)", async () => {
      await session.step(78, "Then \"Show Empty Categories\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Empty Categories", el("box plot viewer"), "true"));
      await session.step(79, "When user zooms into the value axis of box plot viewer", () => zoomValueAxis(page, el("box plot viewer")));
      await session.step(80, "Then box plot viewer should show a narrower value range than before", () => narrowerRange(page, el("box plot viewer")));
      await session.step(81, "When user remembers the value range of box plot viewer", () => rememberRange(page, el("box plot viewer")));
      await session.step(82, "And user sets \"Show Empty Categories\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Empty Categories", el("box plot viewer"), "false"));
      await session.step(83, "Then box plot viewer should show the remembered value range", () => rememberedRange(page, el("box plot viewer")));
      await session.step(84, "When user sets \"Show Empty Categories\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Empty Categories", el("box plot viewer"), "true"));
      await session.step(85, "Then box plot viewer should show the remembered value range", () => rememberedRange(page, el("box plot viewer")));
      await session.step(86, "When user double-clicks on the \"view\" area of box plot viewer", () => doubleClickArea(page, "view", el("box plot viewer")));
      await session.step(87, "Then box plot viewer should show a wider value range than before", () => widerRange(page, el("box plot viewer")));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
