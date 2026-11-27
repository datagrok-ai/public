/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/06-dashboards-and-sharing/08-presentation-mode.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {savedInPresentationMode} from '../../../bindings/howto-dashboards.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, hoverOver, shouldBe, shouldBeSwitchedOn, shouldContainText, switchOn, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {balloonsClosed, browsePanelClosed, browsePanelOpen, closeAllViews, contextPanelClosed, dialogCloses, doubleClickOpensProject, noProjectOnServer, openDataset, simpleModeOff, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, infoBalloonText, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Dashboards", () => {
  const session = feature(test, "features/howto/06-dashboards-and-sharing/08-presentation-mode.feature", import.meta.url);
  test("Open a dashboard in presentation mode", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And simple mode is off", () => simpleModeOff(page));
    await session.step(24, "And no project named \"Demographics overview\" is on the server", () => noProjectOnServer(page, "Demographics overview"));
    await session.step(25, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(26, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]);
    await session.step(30, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","HEIGHT"]]), [["valueColumnName","HEIGHT"]]);
    await session.step(32, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","RACE"]]), [["splitColumnName","RACE"]]);
    await session.step(34, "And user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(35, "And user enters \"Demographics overview\" into Name text input in \"Save project\" dialog", () => enterInto(page, "Demographics overview", el("Name text input in \"Save project\" dialog")));
    await session.step(36, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(37, "And the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(38, "And \"Share Demographics overview\" dialog should be visible", () => shouldBe(page, el("\"Share Demographics overview\" dialog"), "visible"));
    await session.step(39, "And user clicks on CANCEL button in \"Share Demographics overview\" dialog", () => clickOn(page, el("CANCEL button in \"Share Demographics overview\" dialog")));
    await session.step(40, "And the \"Share Demographics overview\" dialog should close", () => dialogCloses(page, "Share Demographics overview"));
    await session.step(41, "And user closes all views", () => closeAllViews(page));
    await session.step(42, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(43, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(44, "And user enters \"Demographics overview\" into gallery search", () => enterInto(page, "Demographics overview", el("gallery search")));
    await session.step(45, "And \"Demographics overview\" gallery card should be visible", () => shouldBe(page, el("\"Demographics overview\" gallery card"), "visible"));
    await session.step(46, "And user double-clicks on \"Demographics overview\" gallery card and the project opens", () => doubleClickOpensProject(page, el("\"Demographics overview\" gallery card")));
    await session.step(47, "And the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(48, "And the notification balloons are closed", () => balloonsClosed(page));
    await session.step(49, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(50, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(51, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(52, "When the video says \"Some people only need to look at the charts\"", () => videoSays(page, "Some people only need to look at the charts"));
    await session.step(53, "And user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(54, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(56, "When user hovers over \"Presentation mode\" input in \"Save project\" dialog", () => hoverOver(page, el("\"Presentation mode\" input in \"Save project\" dialog")), undefined, "Presentation mode: open it with only the charts showing");
    await session.step(57, "Then tooltip should contain text \"Design mode\"", () => shouldContainText(page, el("tooltip"), "Design mode"));
    await session.step(58, "When user switches on \"Presentation mode\" input in \"Save project\" dialog", () => switchOn(page, el("\"Presentation mode\" input in \"Save project\" dialog")));
    await session.step(59, "Then \"Presentation mode\" input in \"Save project\" dialog should be switched on", () => shouldBeSwitchedOn(page, el("\"Presentation mode\" input in \"Save project\" dialog")));
    await session.step(60, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(61, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(62, "And the project named \"Demographics overview\" should be saved to open in presentation mode", () => savedInPresentationMode(page, "Demographics overview"));
    await session.step(63, "Given user closes all views", () => closeAllViews(page));
    await session.step(64, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(65, "When the video says \"Open it again from Browse > Dashboards\"", () => videoSays(page, "Open it again from Browse > Dashboards"));
    await session.step(66, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(67, "And user enters \"Demographics overview\" into gallery search", () => enterInto(page, "Demographics overview", el("gallery search")));
    await session.step(68, "Then \"Demographics overview\" gallery card should be visible", () => shouldBe(page, el("\"Demographics overview\" gallery card"), "visible"));
    await session.step(70, "When user double-clicks on \"Demographics overview\" gallery card and the project opens", () => doubleClickOpensProject(page, el("\"Demographics overview\" gallery card")), undefined, "It opens without the side panels and the status bar");
    await session.step(71, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(72, "And \"back to design mode\" link should be visible", () => shouldBe(page, el("\"back to design mode\" link"), "visible"));
    await session.step(73, "And browse tab should be hidden", () => shouldBe(page, el("browse tab"), "hidden"));
    await session.step(74, "And toolbox tab should be hidden", () => shouldBe(page, el("toolbox tab"), "hidden"));
    await session.step(75, "And status bar should be hidden", () => shouldBe(page, el("status bar"), "hidden"));
    await session.step(76, "And an info balloon containing \"Press F7 to go back to the design mode\" should have been shown", () => infoBalloonText(page, "Press F7 to go back to the design mode"));
    await session.step(78, "When user clicks on \"back to design mode\" link", () => clickOn(page, el("\"back to design mode\" link")), undefined, "back to design mode, or F7, brings the panels back");
    await session.step(79, "Then status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
    await session.step(80, "And browse tab should be visible", () => shouldBe(page, el("browse tab"), "visible"));
    await session.step(81, "And no errors should have been logged", () => noErrors(page));
  });
});
