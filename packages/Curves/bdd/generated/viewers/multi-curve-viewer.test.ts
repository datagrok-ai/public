/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/multi-curve-viewer.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/curves.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, expand, hoverOver, pressKey, shouldBe, shouldContainText, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {clearSelection, selectFirstRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, closeAllViews, contextPanelShowsCurrentCell, openDataset, packageInstalled, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, noErrors, readingIs, readingReads, reportsNoError, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickFromViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Multi curve viewer from the Add Viewer dialog, its curve columns, Use in Trellis", () => {
  const session = feature(test, "features/viewers/multi-curve-viewer.feature", import.meta.url);
  test("The Multi curve viewer from the Add Viewer dialog, its curve columns, Use in Trellis", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And simple mode is off", () => simpleModeOff(page));
    await session.step(22, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(23, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await run.scenario("On a table without curves the MultiCurveViewer card is disabled", async () => {
      await session.step(26, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(27, "Then \"Add viewer\" icon in toolbar should be visible", () => shouldBe(page, el("\"Add viewer\" icon in toolbar"), "visible"));
      await session.step(28, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
      await session.step(29, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
      await session.step(30, "When user types \"curve\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "curve", el("viewer gallery search in \"Add Viewer\" dialog")));
      await session.step(31, "Then first \"MultiCurveViewer\" card in \"Add Viewer\" dialog should be disabled", () => shouldBe(page, el("first \"MultiCurveViewer\" card in \"Add Viewer\" dialog"), "disabled"));
      await session.step(32, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(33, "Then \"Add Viewer\" dialog should be absent", () => shouldBe(page, el("\"Add Viewer\" dialog"), "absent"));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("On the curves table the card adds the viewer with one curve column, and All takes five", async () => {
      await session.step(37, "Given user closes all views", () => closeAllViews(page));
      await session.step(38, "And user opens curves dataset", () => openDataset(page, ds("curves")));
      await session.step(39, "Then the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
      await session.step(40, "And \"Add viewer\" icon in toolbar should be visible", () => shouldBe(page, el("\"Add viewer\" icon in toolbar"), "visible"));
      await session.step(41, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
      await session.step(42, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
      await session.step(43, "When user types \"curve\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "curve", el("viewer gallery search in \"Add Viewer\" dialog")));
      await session.step(44, "Then first \"MultiCurveViewer\" card in \"Add Viewer\" dialog should be enabled", () => shouldBe(page, el("first \"MultiCurveViewer\" card in \"Add Viewer\" dialog"), "enabled"));
      await session.step(45, "When user clicks on first \"MultiCurveViewer\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"MultiCurveViewer\" card in \"Add Viewer\" dialog")));
      await session.step(46, "Then \"Add Viewer\" dialog should be absent", () => shouldBe(page, el("\"Add Viewer\" dialog"), "absent"));
      await session.step(47, "And the open tableview should have 1 MultiCurveViewer viewer", () => viewerCount(page, 1, "MultiCurveViewer"));
      await session.step(48, "And the \"curves\" reading of MultiCurveViewer viewer should be 2", () => readingIs(page, "curves", el("MultiCurveViewer viewer"), 2));
      await session.step(49, "When user picks \"Properties...\" from the viewer menu of MultiCurveViewer viewer", () => pickFromViewerMenu(page, "Properties...", el("MultiCurveViewer viewer")));
      await session.step(50, "Then \"Curves\" property in context panel should contain text \"1 / 5\"", () => shouldContainText(page, el("\"Curves\" property in context panel"), "1 / 5"));
      await session.step(51, "When user clicks on \"...\" button in \"Curves\" property in context panel", () => clickOn(page, el("\"...\" button in \"Curves\" property in context panel")));
      await session.step(52, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(53, "When user clicks on \"All\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"All\" link in \"Select columns...\" dialog")));
      await session.step(54, "Then \"5 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"5 checked\" text in \"Select columns...\" dialog"), "visible"));
      await session.step(55, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(56, "And user picks \"Properties...\" from the viewer menu of MultiCurveViewer viewer", () => pickFromViewerMenu(page, "Properties...", el("MultiCurveViewer viewer")));
      await session.step(57, "Then \"Curves\" property in context panel should contain text \"5 / 5\"", () => shouldContainText(page, el("\"Curves\" property in context panel"), "5 / 5"));
      await session.step(58, "And the \"curves\" reading of MultiCurveViewer viewer should be 7", () => readingIs(page, "curves", el("MultiCurveViewer viewer"), 7));
      await session.step(59, "And the \"points outside the axes\" reading of MultiCurveViewer viewer should be 0", () => readingIs(page, "points outside the axes", el("MultiCurveViewer viewer"), 0));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The selected rows add their curves, Merge Cell Series draws one curve per cell", async () => {
      await session.step(63, "When user selects the first 3 rows", () => selectFirstRows(page, 3));
      await session.step(64, "Then the \"curves\" reading of MultiCurveViewer viewer should be 21", () => readingIs(page, "curves", el("MultiCurveViewer viewer"), 21));
      await session.step(65, "When user expands \"Misc\" category in context panel", () => expand(page, el("\"Misc\" category in context panel")));
      await session.step(66, "And user checks \"Merge Cell Series\" property in context panel", () => check(page, el("\"Merge Cell Series\" property in context panel")));
      await session.step(67, "Then the \"curves\" reading of MultiCurveViewer viewer should be 15", () => readingIs(page, "curves", el("MultiCurveViewer viewer"), 15));
      await session.step(68, "When user unchecks \"Merge Cell Series\" property in context panel", () => uncheck(page, el("\"Merge Cell Series\" property in context panel")));
      await session.step(69, "Then the \"curves\" reading of MultiCurveViewer viewer should be 21", () => readingIs(page, "curves", el("MultiCurveViewer viewer"), 21));
      await session.step(70, "When user unchecks \"Show Selected Rows Curves\" property in context panel", () => uncheck(page, el("\"Show Selected Rows Curves\" property in context panel")));
      await session.step(71, "Then the \"curves\" reading of MultiCurveViewer viewer should be 7", () => readingIs(page, "curves", el("MultiCurveViewer viewer"), 7));
      await session.step(72, "When user checks \"Show Selected Rows Curves\" property in context panel", () => check(page, el("\"Show Selected Rows Curves\" property in context panel")));
      await session.step(73, "Then the \"curves\" reading of MultiCurveViewer viewer should be 21", () => readingIs(page, "curves", el("MultiCurveViewer viewer"), 21));
      await session.step(74, "When user clears the row selection", () => clearSelection(page));
      await session.step(75, "Then the \"curves\" reading of MultiCurveViewer viewer should be 7", () => readingIs(page, "curves", el("MultiCurveViewer viewer"), 7));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A closed viewer leaves no error when rows are clicked afterwards", async () => {
      await session.step(79, "When user hovers over MultiCurveViewer viewer", () => hoverOver(page, el("MultiCurveViewer viewer")));
      await session.step(80, "And user clicks on close icon of MultiCurveViewer viewer", () => clickOn(page, el("close icon of MultiCurveViewer viewer")));
      await session.step(81, "Then the open tableview should have 0 MultiCurveViewer viewers", () => viewerCount(page, 0, "MultiCurveViewer"));
      await session.step(82, "When user selects the first 3 rows", () => selectFirstRows(page, 3));
      await session.step(83, "And user makes row 4 current", () => makeRowCurrent(page, 4));
      await session.step(84, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Use in Trellis puts the viewer into a trellis plot as its inner viewer", async () => {
      await session.step(87, "Given user adds MultiCurveViewer viewer", () => addViewer(page, "MultiCurveViewer"));
      await session.step(88, "When user picks \"General | Use in Trellis\" from the viewer menu of MultiCurveViewer viewer", () => pickFromViewerMenu(page, "General | Use in Trellis", el("MultiCurveViewer viewer")));
      await session.step(89, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(90, "And the \"inner viewer type\" reading of trellis plot viewer should be \"MultiCurveViewer\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "MultiCurveViewer"));
      await session.step(91, "And trellis plot viewer should report no error", () => reportsNoError(page, el("trellis plot viewer")));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
