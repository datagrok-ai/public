/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/sunburst/sunburst-inherit-from-grid.feature
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
import {check, clickOn, isExpanded, shouldBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {colorCategorical, colorCodedCategorically} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaColor, areaNotColor, noErrors, pickFromAreaContextMenu, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Sunburst takes its segment colors from the grid's color coding with Inherit From Grid", () => {
  const session = feature(test, "features/viewers/sunburst/sunburst-inherit-from-grid.feature", import.meta.url);
  test("Segments follow the grid's categorical colors, and their own palette with the option off", {tag: ["@viewers", "@realizes:charts.viewer.sunburst"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(16, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(17, "When user clicks on first \"Sunburst\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Sunburst\" card in \"Add Viewer\" dialog")));
    await session.step(18, "Then sunburst viewer should be visible", () => shouldBe(page, el("sunburst viewer"), "visible"));
    await session.step(19, "When user picks \"Properties...\" from the context menu of the \"view\" area of sunburst viewer", () => pickFromAreaContextMenu(page, "Properties...", "view", el("sunburst viewer")));
    await session.step(20, "And user clicks on \"...\" button in \"Hierarchy\" property in context panel", () => clickOn(page, el("\"...\" button in \"Hierarchy\" property in context panel")));
    await session.step(21, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(22, "When user clicks on \"None\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"None\" link in \"Select columns...\" dialog")));
    await session.step(23, "Then \"0 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"0 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(24, "When user toggles the \"SEX\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "SEX", el("\"Select columns...\" dialog")));
    await session.step(25, "Then \"1 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"1 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(26, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(27, "Then \"Select columns...\" dialog should be absent", () => shouldBe(page, el("\"Select columns...\" dialog"), "absent"));
    await session.step(28, "And the \"hierarchy columns\" reading of sunburst viewer should be \"SEX\"", () => readingReads(page, "hierarchy columns", el("sunburst viewer"), "SEX"));
    await session.step(29, "And the \"segments\" reading of sunburst viewer should be 2", () => readingIs(page, "segments", el("sunburst viewer"), 2));
    await session.step(32, "Then the \"segment F\" area of sunburst viewer should contain the color \"#1F77B4\"", () => areaColor(page, "segment F", el("sunburst viewer"), "#1F77B4"));
    await session.step(33, "And the \"segment F\" area of sunburst viewer should not contain the color \"#FF0000\"", () => areaNotColor(page, "segment F", el("sunburst viewer"), "#FF0000"));
    await session.step(34, "When user colors \"SEX\" column categorically:", () => colorCategorical(page, "SEX", [["F","#FF0000"],["M","#0000FF"]]), [["F","#FF0000"],["M","#0000FF"]]);
    await session.step(37, "Then \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
    await session.step(38, "And the \"segment F\" area of sunburst viewer should contain the color \"#FF0000\"", () => areaColor(page, "segment F", el("sunburst viewer"), "#FF0000"));
    await session.step(39, "And the \"segment M\" area of sunburst viewer should contain the color \"#0000FF\"", () => areaColor(page, "segment M", el("sunburst viewer"), "#0000FF"));
    await session.step(40, "And the \"segment F\" area of sunburst viewer should not contain the color \"#1F77B4\"", () => areaNotColor(page, "segment F", el("sunburst viewer"), "#1F77B4"));
    await session.step(41, "Given \"Color\" category in context panel is expanded", () => isExpanded(page, el("\"Color\" category in context panel")));
    await session.step(42, "When user unchecks \"Inherit From Grid\" property in context panel", () => uncheck(page, el("\"Inherit From Grid\" property in context panel")));
    await session.step(43, "Then the \"segment F\" area of sunburst viewer should contain the color \"#1F77B4\"", () => areaColor(page, "segment F", el("sunburst viewer"), "#1F77B4"));
    await session.step(44, "And the \"segment F\" area of sunburst viewer should not contain the color \"#FF0000\"", () => areaNotColor(page, "segment F", el("sunburst viewer"), "#FF0000"));
    await session.step(45, "And the \"segment M\" area of sunburst viewer should not contain the color \"#0000FF\"", () => areaNotColor(page, "segment M", el("sunburst viewer"), "#0000FF"));
    await session.step(46, "When user checks \"Inherit From Grid\" property in context panel", () => check(page, el("\"Inherit From Grid\" property in context panel")));
    await session.step(47, "Then the \"segment F\" area of sunburst viewer should contain the color \"#FF0000\"", () => areaColor(page, "segment F", el("sunburst viewer"), "#FF0000"));
    await session.step(48, "And the \"segment M\" area of sunburst viewer should contain the color \"#0000FF\"", () => areaColor(page, "segment M", el("sunburst viewer"), "#0000FF"));
    await session.step(49, "And no errors should have been logged", () => noErrors(page));
  });
});
