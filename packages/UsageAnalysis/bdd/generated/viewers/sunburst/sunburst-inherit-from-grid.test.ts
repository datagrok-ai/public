/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/sunburst/sunburst-inherit-from-grid.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.sunburst, public#2292]
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
import {check, clickOn, isExpanded, selectIn, shouldBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {colorCategorical, colorCodedCategorically} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaColor, areaNotColor, boundTable, noErrors, pickFromAreaContextMenu, readingIs, readingNotAsRemembered, readingReads, rememberReading, reportsNoError} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Sunburst takes its segment colors from the grid's color coding with Inherit From Grid, and follows its Table property", () => {
  const session = feature(test, "features/viewers/sunburst/sunburst-inherit-from-grid.feature", import.meta.url);
  test("Sunburst takes its segment colors from the grid's color coding with Inherit From Grid, and follows its Table property", {tag: ["@journey", "@viewers", "@realizes:charts.viewer.sunburst", "@realizes:public#2292"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(19, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(20, "When user clicks on first \"Sunburst\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Sunburst\" card in \"Add Viewer\" dialog")));
    await session.step(21, "Then sunburst viewer should be visible", () => shouldBe(page, el("sunburst viewer"), "visible"));
    await session.step(22, "When user picks \"Properties...\" from the context menu of the \"view\" area of sunburst viewer", () => pickFromAreaContextMenu(page, "Properties...", "view", el("sunburst viewer")));
    await session.step(23, "And user clicks on \"...\" button in \"Hierarchy\" property in context panel", () => clickOn(page, el("\"...\" button in \"Hierarchy\" property in context panel")));
    await session.step(24, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(25, "When user clicks on \"None\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"None\" link in \"Select columns...\" dialog")));
    await session.step(26, "Then \"0 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"0 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(27, "When user toggles the \"SEX\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "SEX", el("\"Select columns...\" dialog")));
    await session.step(28, "Then \"1 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"1 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(29, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(30, "Then \"Select columns...\" dialog should be absent", () => shouldBe(page, el("\"Select columns...\" dialog"), "absent"));
    await session.step(31, "And the \"hierarchy columns\" reading of sunburst viewer should be \"SEX\"", () => readingReads(page, "hierarchy columns", el("sunburst viewer"), "SEX"));
    await session.step(32, "And the \"segments\" reading of sunburst viewer should be 2", () => readingIs(page, "segments", el("sunburst viewer"), 2));
    await run.scenario("Segments follow the grid's categorical colors, and their own palette with the option off", async () => {
      await session.step(35, "Then the \"segment F\" area of sunburst viewer should contain the color \"#1F77B4\"", () => areaColor(page, "segment F", el("sunburst viewer"), "#1F77B4"));
      await session.step(36, "And the \"segment F\" area of sunburst viewer should not contain the color \"#FF0000\"", () => areaNotColor(page, "segment F", el("sunburst viewer"), "#FF0000"));
      await session.step(37, "When user colors \"SEX\" column categorically:", () => colorCategorical(page, "SEX", [["F","#FF0000"],["M","#0000FF"]]), [["F","#FF0000"],["M","#0000FF"]]);
      await session.step(40, "Then \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(41, "And the \"segment F\" area of sunburst viewer should contain the color \"#FF0000\"", () => areaColor(page, "segment F", el("sunburst viewer"), "#FF0000"));
      await session.step(42, "And the \"segment M\" area of sunburst viewer should contain the color \"#0000FF\"", () => areaColor(page, "segment M", el("sunburst viewer"), "#0000FF"));
      await session.step(43, "And the \"segment F\" area of sunburst viewer should not contain the color \"#1F77B4\"", () => areaNotColor(page, "segment F", el("sunburst viewer"), "#1F77B4"));
      await session.step(44, "Given \"Color\" category in context panel is expanded", () => isExpanded(page, el("\"Color\" category in context panel")));
      await session.step(45, "When user unchecks \"Inherit From Grid\" property in context panel", () => uncheck(page, el("\"Inherit From Grid\" property in context panel")));
      await session.step(46, "Then the \"segment F\" area of sunburst viewer should contain the color \"#1F77B4\"", () => areaColor(page, "segment F", el("sunburst viewer"), "#1F77B4"));
      await session.step(47, "And the \"segment F\" area of sunburst viewer should not contain the color \"#FF0000\"", () => areaNotColor(page, "segment F", el("sunburst viewer"), "#FF0000"));
      await session.step(48, "And the \"segment M\" area of sunburst viewer should contain the color \"#FFBB78\"", () => areaColor(page, "segment M", el("sunburst viewer"), "#FFBB78"));
      await session.step(49, "And the \"segment M\" area of sunburst viewer should not contain the color \"#0000FF\"", () => areaNotColor(page, "segment M", el("sunburst viewer"), "#0000FF"));
      await session.step(50, "When user checks \"Inherit From Grid\" property in context panel", () => check(page, el("\"Inherit From Grid\" property in context panel")));
      await session.step(51, "Then the \"segment F\" area of sunburst viewer should contain the color \"#FF0000\"", () => areaColor(page, "segment F", el("sunburst viewer"), "#FF0000"));
      await session.step(52, "And the \"segment M\" area of sunburst viewer should contain the color \"#0000FF\"", () => areaColor(page, "segment M", el("sunburst viewer"), "#0000FF"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Table property moves the sunburst to another table and back (public#2292)", async () => {
      await session.step(57, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
      await session.step(58, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(59, "When user remembers the \"segment names\" reading of sunburst viewer", () => rememberReading(page, "segment names", el("sunburst viewer")));
      await session.step(60, "And user clicks on settings icon of sunburst viewer", () => clickOn(page, el("settings icon of sunburst viewer")));
      await session.step(61, "And user selects \"spgi-100\" in \"Table\" property in context panel", () => selectIn(page, "spgi-100", el("\"Table\" property in context panel")));
      await session.step(62, "Then sunburst viewer should be bound to table \"spgi-100\"", () => boundTable(page, el("sunburst viewer"), "spgi-100"));
      await session.step(63, "And the \"segment names\" reading of sunburst viewer should not be as remembered", () => readingNotAsRemembered(page, "segment names", el("sunburst viewer")));
      await session.step(64, "And sunburst viewer should report no error", () => reportsNoError(page, el("sunburst viewer")));
      await session.step(65, "When user selects \"demog-1000\" in \"Table\" property in context panel", () => selectIn(page, "demog-1000", el("\"Table\" property in context panel")));
      await session.step(66, "Then sunburst viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("sunburst viewer"), "demog-1000"));
      await session.step(67, "And the \"hierarchy columns\" reading of sunburst viewer should be \"SEX, CONTROL, RACE\"", () => readingReads(page, "hierarchy columns", el("sunburst viewer"), "SEX, CONTROL, RACE"));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
