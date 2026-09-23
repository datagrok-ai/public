/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scripts/scripts-layout.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.scripts]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {cleanLayouts, saveScript} from '../../bindings/scripts.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {entityHasLayout, scriptOnServer, scriptsView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, painted, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The layout of a script's result", () => {
  const session = feature(test, "features/scripts/scripts-layout.feature", import.meta.url);
  test("The layout of a script's result", {tag: ["@journey", "@serial", "@realizes:views.scripts"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And a script \"BddScriptLayout{time}\" is on the server:", () => scriptOnServer(page, session.text("BddScriptLayout{time}"), "//language: javascript\n//output: dataframe df\ndf = await grok.data.getDemoTable('cars.csv');"));
    await session.step(25, "And the layouts saved for the script are deleted at the end", () => cleanLayouts(page));
    await session.step(26, "And user opens the Scripts view", () => scriptsView(page));
    await run.scenario("The Layout tab asks for a run first", async () => {
      await session.step(29, "When user types \"BddScriptLayout{time}\" into gallery search", () => typeInto(page, session.text("BddScriptLayout{time}"), el("gallery search")));
      await session.step(30, "And user double-clicks on \"BddScriptLayout{time}\" link in gallery", () => doubleClickOn(page, el(session.text("\"BddScriptLayout{time}\" link in gallery"))));
      await session.step(31, "Then the \"BddScriptLayout{time}\" view should be current", () => viewIsCurrent(page, session.text("BddScriptLayout{time}")));
      await session.step(32, "When user clicks on Layout tab", () => clickOn(page, el("Layout tab")));
      await session.step(33, "Then \"Run script to get data and edit layout.\" text should be visible", () => shouldBe(page, el("\"Run script to get data and edit layout.\" text"), "visible"));
      await session.step(34, "And grid viewer should be absent", () => shouldBe(page, el("grid viewer"), "absent"));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Running the script fills the layout with its result", async () => {
      await session.step(38, "When user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")));
      await session.step(39, "Then grid viewer should be visible", () => shouldBe(page, el("grid viewer"), "visible"));
      await session.step(40, "And the \"rows\" reading of grid should be 30", () => readingIs(page, "rows", el("grid"), 30));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
      await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A viewer added from the toolbox joins the layout", async () => {
      await session.step(45, "When user clicks on \"bar chart\" icon in toolbox", () => clickOn(page, el("\"bar chart\" icon in toolbox")));
      await session.step(46, "Then bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
      await session.step(47, "And bar chart viewer should be painted", () => painted(page, el("bar chart viewer")));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
      await session.step(49, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Save stores the layout with the script", async () => {
      await session.step(52, "When user saves the script", () => saveScript(page));
      await session.step(53, "Then no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(54, "And the script \"BddScriptLayout{time}\" on the server should have a layout", () => entityHasLayout(page, "script", session.text("BddScriptLayout{time}")));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
