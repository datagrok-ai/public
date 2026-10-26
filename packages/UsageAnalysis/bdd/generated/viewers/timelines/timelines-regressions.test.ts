/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/timelines/timelines-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.timelines, public#617, public#1973]
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
import {clickOn, doubleClickOn, isExpanded, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCalculated, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelOpen, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {legendLists, noErrors, painted, propertyShouldBe, readingIs, repainted, reportsNoError, setProperty, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Timelines regressions of closed issues", () => {
  const session = feature(test, "features/viewers/timelines/timelines-regressions.feature", import.meta.url);
  test("Timelines regressions of closed issues", {tag: ["@journey", "@viewers", "@realizes:charts.viewer.timelines", "@realizes:public#617", "@realizes:public#1973"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(15, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(16, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(17, "And Files---App-Data---Charts tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Charts tree node inside browse tree")));
    await session.step(18, "When user double-clicks on Files---App-Data---Charts---ae.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Charts---ae.csv tree node inside browse tree")));
    await session.step(19, "Then the \"ae\" view should be current", () => viewIsCurrent(page, "ae"));
    await session.step(20, "And the table should have 143 rows", () => rowCount(page, 143));
    await session.step(21, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(22, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(23, "When user clicks on first \"Timelines\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Timelines\" card in \"Add Viewer\" dialog")));
    await session.step(24, "Then timelines viewer should be visible", () => shouldBe(page, el("timelines viewer"), "visible"));
    await session.step(25, "And timelines viewer should be painted", () => painted(page, el("timelines viewer")));
    await session.step(26, "And the \"intervals\" reading of timelines viewer should be 143", () => readingIs(page, "intervals", el("timelines viewer"), 143));
    await session.step(27, "When user clicks on grid", () => clickOn(page, el("grid")));
    await session.step(28, "And user clicks on settings icon of timelines viewer", () => clickOn(page, el("settings icon of timelines viewer")));
    await run.scenario("Clearing Color keeps every interval and throws nothing (public#617)", async () => {
      await session.step(32, "Given \"Color\" category in context panel is expanded", () => isExpanded(page, el("\"Color\" category in context panel")));
      await session.step(33, "When user selects \"AESOC\" in \"Color\" property in context panel", () => selectIn(page, "AESOC", el("\"Color\" property in context panel")));
      await session.step(34, "Then the legend of timelines viewer should list 15 items", () => legendLists(page, el("timelines viewer"), 15));
      await session.step(35, "When user takes a snapshot of timelines viewer", () => takeSnapshot(page, el("timelines viewer")));
      await session.step(36, "And user sets \"colorColumnName\" property of timelines viewer to \"\"", () => setProperty(page, "colorColumnName", el("timelines viewer"), ""));
      await session.step(37, "Then timelines viewer should have repainted", () => repainted(page, el("timelines viewer")));
      await session.step(38, "And legend of timelines viewer should be hidden", () => shouldBe(page, el("legend of timelines viewer"), "hidden"));
      await session.step(39, "And the \"intervals\" reading of timelines viewer should be 143", () => readingIs(page, "intervals", el("timelines viewer"), 143));
      await session.step(40, "And timelines viewer should report no error", () => reportsNoError(page, el("timelines viewer")));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Long Split By categories still draw their lanes (public#1973)", async () => {
      await session.step(45, "Given \"Data\" category in context panel is expanded", () => isExpanded(page, el("\"Data\" category in context panel")));
      await session.step(46, "Then \"Split By\" property of timelines viewer should be \"USUBJID\"", () => propertyShouldBe(page, "Split By", el("timelines viewer"), "USUBJID"));
      await session.step(47, "And the \"lanes\" reading of timelines viewer should be 71", () => readingIs(page, "lanes", el("timelines viewer"), 71));
      await session.step(48, "When user adds a calculated column \"LongSubj\" with formula \"${USUBJID} + ' is a deliberately very long subject label that runs past the axis'\"", () => addCalculated(page, "LongSubj", "${USUBJID} + ' is a deliberately very long subject label that runs past the axis'"));
      await session.step(49, "And user selects \"LongSubj\" in \"Split By\" property in context panel", () => selectIn(page, "LongSubj", el("\"Split By\" property in context panel")));
      await session.step(50, "Then \"Split By\" property of timelines viewer should be \"LongSubj\"", () => propertyShouldBe(page, "Split By", el("timelines viewer"), "LongSubj"));
      await session.step(51, "And the \"lanes\" reading of timelines viewer should be 71", () => readingIs(page, "lanes", el("timelines viewer"), 71));
      await session.step(52, "And the \"intervals\" reading of timelines viewer should be 143", () => readingIs(page, "intervals", el("timelines viewer"), 143));
      await session.step(53, "And timelines viewer should report no error", () => reportsNoError(page, el("timelines viewer")));
      await session.step(54, "And timelines viewer should be painted", () => painted(page, el("timelines viewer")));
      await session.step(55, "When user selects \"USUBJID\" in \"Split By\" property in context panel", () => selectIn(page, "USUBJID", el("\"Split By\" property in context panel")));
      await session.step(56, "Then \"Split By\" property of timelines viewer should be \"USUBJID\"", () => propertyShouldBe(page, "Split By", el("timelines viewer"), "USUBJID"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
