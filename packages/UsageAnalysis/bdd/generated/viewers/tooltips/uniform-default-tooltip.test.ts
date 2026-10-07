/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tooltips/uniform-default-tooltip.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tooltips]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/grid.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hoverArea, noErrors, pointerAway, setProperties, tooltipColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Viewers show the same default tooltip", () => {
  const session = feature(test, "features/viewers/tooltips/uniform-default-tooltip.feature", import.meta.url);
  test("Viewers show the same default tooltip", {tag: ["@journey", "@viewers", "@realizes:viewers.tooltips"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(24, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(25, "And user adds a box plot viewer", () => addViewer(page, "box plot"));
    await session.step(26, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await session.step(28, "And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await session.step(30, "And user sets properties of grid:", () => setProperties(page, el("grid"), [["Show Tooltip","inherit from table"],["Show Column Names","Always"],["Show Visible Columns In Tooltip","true"]]), [["Show Tooltip","inherit from table"],["Show Column Names","Always"],["Show Visible Columns In Tooltip","true"]]);
    await run.scenario("The scatter plot, the box plot and the grid list the same columns", async () => {
      await session.step(36, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(37, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(38, "And the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(39, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
      await session.step(40, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(41, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(42, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(43, "And the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(44, "When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
      await session.step(45, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(46, "When user hovers over the \"cell 11 of AGE\" area of grid", () => hoverArea(page, "cell 11 of AGE", el("grid")));
      await session.step(47, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(48, "And the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(49, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(50, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
