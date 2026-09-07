/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/box-plot/box-plot-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {addCalculated, filterBetween, filterIsExactly, filterPassesAll, filterPassesFewer, removeColumn, resetFilter, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, hasNoArea, narrowerRange, noErrors, propertyShouldBe, repainted, sameRange, setProperty, widerRange} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Box plot filter semantics", () => {
  const session = feature(test, "features/viewers/box-plot/box-plot-filter.feature", import.meta.url);
  test("Box plot filter semantics", {tag: ["@journey", "@viewers", "@realizes:viewers.box-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(11, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","Average Mass"],["Category 1","Series"]]));
    await session.step(14, "Then the table should have 100 rows", () => rowCount(page, 100));
    await run.scenario("The value range follows the filter", async () => {
      await session.step(17, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(18, "And \"Zoom Values By Filter\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Zoom Values By Filter", el("box plot viewer"), "true"));
      await session.step(19, "When user filters rows where \"Average Mass\" is between 300 and 400", () => filterBetween(page, "Average Mass", 300, 400));
      await session.step(20, "Then fewer than 100 rows should pass the filter", () => filterPassesFewer(page, 100));
      await session.step(21, "And box plot viewer should show a narrower value range than before", () => narrowerRange(page, el("box plot viewer")));
      await session.step(22, "When user resets the filter", () => resetFilter(page));
      await session.step(23, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(24, "And box plot viewer should show a wider value range than before", () => widerRange(page, el("box plot viewer")));
    });
    await run.scenario("Zoom Values By Filter off keeps the range", async () => {
      await session.step(27, "When user sets \"Zoom Values By Filter\" property of box plot viewer to \"false\"", () => setProperty(page, "Zoom Values By Filter", el("box plot viewer"), "false"));
      await session.step(28, "And user filters rows where \"Average Mass\" is between 300 and 400", () => filterBetween(page, "Average Mass", 300, 400));
      await session.step(29, "Then fewer than 100 rows should pass the filter", () => filterPassesFewer(page, 100));
      await session.step(30, "And box plot viewer should show the same value range as before", () => sameRange(page, el("box plot viewer")));
      await session.step(31, "When user sets \"Zoom Values By Filter\" property of box plot viewer to \"true\"", () => setProperty(page, "Zoom Values By Filter", el("box plot viewer"), "true"));
      await session.step(32, "Then box plot viewer should show a narrower value range than before", () => narrowerRange(page, el("box plot viewer")));
      await session.step(33, "When user resets the filter", () => resetFilter(page));
    });
    await run.scenario("The viewer's own filter leaves the table's alone", async () => {
      await session.step(36, "When user filters rows where \"Average Mass\" is between 300 and 400", () => filterBetween(page, "Average Mass", 300, 400));
      await session.step(37, "And user sets \"Filter\" property of box plot viewer to \"${Average Mass} > 350\"", () => setProperty(page, "Filter", el("box plot viewer"), "${Average Mass} > 350"));
      await session.step(38, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(39, "And the filter should pass exactly the rows where \"Average Mass\" is between 300 and 400", () => filterIsExactly(page, "Average Mass", 300, 400));
      await session.step(40, "When user sets \"Filter\" property of box plot viewer to \"\"", () => setProperty(page, "Filter", el("box plot viewer"), ""));
      await session.step(41, "And user resets the filter", () => resetFilter(page));
    });
    await run.scenario("Show Empty Categories drops and restores an empty-valued category", async () => {
      await session.step(44, "When user adds a calculated column \"AverageMassFixture\" with formula \"if(${Series} == \\\"Triazoles\\\", null, ${Average Mass})\"", () => addCalculated(page, "AverageMassFixture", "if(${Series} == \"Triazoles\", null, ${Average Mass})"));
      await session.step(45, "And user sets \"Value\" property of box plot viewer to \"AverageMassFixture\"", () => setProperty(page, "Value", el("box plot viewer"), "AverageMassFixture"));
      await session.step(46, "And user sets \"Show Empty Categories\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Empty Categories", el("box plot viewer"), "true"));
      await session.step(47, "Then box plot viewer should have a \"category Triazoles\" area", () => hasArea(page, el("box plot viewer"), "category Triazoles"));
      await session.step(48, "When user sets \"Show Empty Categories\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Empty Categories", el("box plot viewer"), "false"));
      await session.step(49, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(50, "And box plot viewer should not have a \"category Triazoles\" area", () => hasNoArea(page, el("box plot viewer"), "category Triazoles"));
      await session.step(51, "When user sets \"Show Empty Categories\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Empty Categories", el("box plot viewer"), "true"));
      await session.step(52, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(53, "And box plot viewer should have a \"category Triazoles\" area", () => hasArea(page, el("box plot viewer"), "category Triazoles"));
      await session.step(54, "When user sets \"Value\" property of box plot viewer to \"Average Mass\"", () => setProperty(page, "Value", el("box plot viewer"), "Average Mass"));
      await session.step(55, "And user removes \"AverageMassFixture\" column", () => removeColumn(page, "AverageMassFixture"));
    });
    await run.scenario("A coloring under a filter keeps its color scale", async () => {
      await session.step(58, "When user filters rows where \"Average Mass\" is between 300 and 400", () => filterBetween(page, "Average Mass", 300, 400));
      await session.step(59, "And user sets \"Marker Color Column\" property of box plot viewer to \"TPSA\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "TPSA"));
      await session.step(60, "Then \"Marker Color Column\" property of box plot viewer should be \"TPSA\"", () => propertyShouldBe(page, "Marker Color Column", el("box plot viewer"), "TPSA"));
      await session.step(61, "And box plot viewer should have a \"color scale\" area", () => hasArea(page, el("box plot viewer"), "color scale"));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
      await session.step(63, "When user resets the filter", () => resetFilter(page));
      await session.step(64, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(65, "And box plot viewer should have a \"color scale\" area", () => hasArea(page, el("box plot viewer"), "color scale"));
      await session.step(66, "When user sets \"Marker Color Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), ""));
    });
    run.finish();
  });
});
