/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/curves-in-trellis.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/curves.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, expand, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, closeAllViews, openDataset, packageInstalled, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, clickArea, hoverArea, loadLayout, noErrors, readingIs, readingReads, reportsNoError, saveLayoutToServer, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {innerPropertyShouldBe, noSuchReading, pickInnerViewer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Multi curve viewer inside a Trellis plot", () => {
  const session = feature(test, "features/viewers/curves-in-trellis.feature", import.meta.url);
  test("The Multi curve viewer inside a Trellis plot", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And simple mode is off", () => simpleModeOff(page));
    await session.step(21, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(22, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await run.scenario("A Multi curve viewer is created in the session before the trellis", async () => {
      await session.step(25, "Given user opens curves dataset", () => openDataset(page, ds("curves")));
      await session.step(26, "And user adds MultiCurveViewer viewer", () => addViewer(page, "MultiCurveViewer"));
      await session.step(27, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Curves becomes the inner viewer of a trellis split by molecule", async () => {
      await session.step(30, "Given user opens curves dataset", () => openDataset(page, ds("curves")));
      await session.step(31, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","smiles"]]), [["X Column Names","smiles"]]);
      await session.step(33, "When user picks \"Curves\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Curves", el("trellis plot viewer")));
      await session.step(34, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"MultiCurveViewer\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "MultiCurveViewer"));
      await session.step(35, "And trellis plot viewer should report no error", () => reportsNoError(page, el("trellis plot viewer")));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The full-screen view of a cell shows its row and closes with one click on X", async () => {
      await session.step(39, "When user hovers over the \"cell body 1,1\" area of trellis plot viewer", () => hoverArea(page, "cell body 1,1", el("trellis plot viewer")));
      await session.step(40, "And user clicks on the \"full screen icon\" area of trellis plot viewer", () => clickArea(page, "full screen icon", el("trellis plot viewer")));
      await session.step(41, "Then the \"full screen rows\" reading of trellis plot viewer should be 1", () => readingIs(page, "full screen rows", el("trellis plot viewer"), 1));
      await session.step(42, "When user clicks on close icon in dialog", () => clickOn(page, el("close icon in dialog")));
      await session.step(43, "Then trellis plot viewer should not report a \"full screen rows\" reading", () => noSuchReading(page, el("trellis plot viewer"), "full screen rows"));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The MultiCurveViewer tab lists the inner properties, Log Y reaches the trellis, Curves keeps the tab", async () => {
      await session.step(47, "When user clicks on settings icon of trellis plot viewer", () => clickOn(page, el("settings icon of trellis plot viewer")));
      await session.step(48, "And user clicks on \"MultiCurveViewer\" tab in context panel", () => clickOn(page, el("\"MultiCurveViewer\" tab in context panel")));
      await session.step(49, "And user expands \"Misc\" category in context panel", () => expand(page, el("\"Misc\" category in context panel")));
      await session.step(50, "Then \"Merge Column Series\" property in context panel should be visible", () => shouldBe(page, el("\"Merge Column Series\" property in context panel"), "visible"));
      await session.step(51, "When user checks \"Log Y\" property in context panel", () => check(page, el("\"Log Y\" property in context panel")));
      await session.step(52, "Then \"logY\" inner property of trellis plot viewer should be \"true\"", () => innerPropertyShouldBe(page, "logY", el("trellis plot viewer"), "true"));
      await session.step(53, "When user clicks on \"...\" button in \"Curves\" property in context panel", () => clickOn(page, el("\"...\" button in \"Curves\" property in context panel")));
      await session.step(54, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(55, "When user clicks on \"All\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"All\" link in \"Select columns...\" dialog")));
      await session.step(56, "Then \"5 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"5 checked\" text in \"Select columns...\" dialog"), "visible"));
      await session.step(57, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(58, "Then \"Curves\" property in context panel should contain text \"5 / 5\"", () => shouldContainText(page, el("\"Curves\" property in context panel"), "5 / 5"));
      await session.step(59, "And \"Merge Column Series\" property in context panel should be visible", () => shouldBe(page, el("\"Merge Column Series\" property in context panel"), "visible"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A layout saved and applied again keeps Log Y of the inner viewer", async () => {
      await session.step(63, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(64, "And user loads the saved layout", () => loadLayout(page));
      await session.step(65, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(66, "And \"logY\" inner property of trellis plot viewer should be \"true\"", () => innerPropertyShouldBe(page, "logY", el("trellis plot viewer"), "true"));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("On a table without curves every cell says it needs a curve column", async () => {
      await session.step(70, "Given user closes all views", () => closeAllViews(page));
      await session.step(71, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(72, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"]]), [["X Column Names","SEX"]]);
      await session.step(74, "When user picks \"Curves\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Curves", el("trellis plot viewer")));
      await session.step(75, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"MultiCurveViewer\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "MultiCurveViewer"));
      await session.step(76, "When user hovers over the \"cell body F | High\" area of trellis plot viewer", () => hoverArea(page, "cell body F | High", el("trellis plot viewer")));
      await session.step(77, "And user hovers over the \"cell body M | Low\" area of trellis plot viewer", () => hoverArea(page, "cell body M | Low", el("trellis plot viewer")));
      await session.step(78, "Then \"The MultiCurveViewer viewer requires a minimum of 1 curves column.\" text in trellis plot viewer should be visible", () => shouldBe(page, el("\"The MultiCurveViewer viewer requires a minimum of 1 curves column.\" text in trellis plot viewer"), "visible"));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Switching the inner viewer of a trellis on a table without curves logs no error", async () => {
      await session.step(82, "When user picks \"Scatter plot\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Scatter plot", el("trellis plot viewer")));
      await session.step(83, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Scatter plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Scatter plot"));
      await session.step(84, "When user picks \"Bar chart\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Bar chart", el("trellis plot viewer")));
      await session.step(85, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Bar chart\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Bar chart"));
      await session.step(86, "When user picks \"Box plot\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Box plot", el("trellis plot viewer")));
      await session.step(87, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Box plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Box plot"));
      await session.step(88, "When user picks \"Curves\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Curves", el("trellis plot viewer")));
      await session.step(89, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"MultiCurveViewer\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "MultiCurveViewer"));
      await session.step(90, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
