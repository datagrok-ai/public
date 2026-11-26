/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/16-curves/01-curves-in-a-table.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
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
import {hoverCurvePoint, noPointLeftOut} from '../../../bindings/howto-curves.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, expand, isExpanded, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelOpen, contextPanelClosed, packageInstalled, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Curves", () => {
  const session = feature(test, "features/howto/16-curves/01-curves-in-a-table.feature", import.meta.url);
  test("Dose-response curves in a table", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(19, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(20, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(25, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(27, "When user expands Files---App-Data---Curves tree node inside browse tree", () => expand(page, el("Files---App-Data---Curves tree node inside browse tree")), undefined, "Sample curves: Files > App Data > Curves > demo");
    await session.step(28, "And user expands Files---App-Data---Curves---demo tree node inside browse tree", () => expand(page, el("Files---App-Data---Curves---demo tree node inside browse tree")));
    await session.step(29, "And user double-clicks on Files---App-Data---Curves---demo---Assay_Curves.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Curves---demo---Assay_Curves.csv tree node inside browse tree")));
    await session.step(30, "Then the \"Assay_Curves\" view should be current", () => viewIsCurrent(page, "Assay_Curves"));
    await session.step(31, "And the table should have 43 rows", () => rowCount(page, 43));
    await session.step(32, "And \"Fitted Curve\" column should have semantic type \"fit\"", () => columnSemType(page, "Fitted Curve", "fit"));
    await session.step(33, "When the video says \"Each cell draws a dose-response curve: the points and the fit\"", () => videoSays(page, "Each cell draws a dose-response curve: the points and the fit"));
    await session.step(35, "And user hovers over point 12 of the curve in the \"cell 1 of Fitted Curve\" area of grid", () => hoverCurvePoint(page, 12, "cell 1 of Fitted Curve", el("grid")), undefined, "Point at a measurement: its concentration and response");
    await session.step(37, "And user clicks on the \"cell 1 of Fitted Curve\" area of grid", () => clickArea(page, "cell 1 of Fitted Curve", el("grid")), undefined, "Click a curve, then open the context panel");
    await session.step(38, "Then no point of the curve in row 1 of the \"Fitted Curve\" column should be left out", () => noPointLeftOut(page, 1, "Fitted Curve"));
    await session.step(39, "When user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")));
    await session.step(40, "Then \"Chart\" pane in context panel should be visible", () => shouldBe(page, el("\"Chart\" pane in context panel"), "visible"));
    await session.step(42, "When user expands Chart accordion header in context panel", () => expand(page, el("Chart accordion header in context panel")), undefined, "Chart: the same curve, larger, with its axes labelled");
    await session.step(44, "And user expands Fit accordion header in context panel", () => expand(page, el("Fit accordion header in context panel")), undefined, "Fit: IC50, Hill slope, R² and more, fitted in your browser");
    await session.step(45, "Then \"IC50\" input in context panel should be visible", () => shouldBe(page, el("\"IC50\" input in context panel"), "visible"));
    await session.step(47, "When user double-clicks on the \"cell 1 of Fitted Curve\" area of grid", () => doubleClickArea(page, "cell 1 of Fitted Curve", el("grid")), undefined, "Double-click a curve to open it in a window of its own");
    await session.step(48, "Then \"Edit chart\" dialog should be visible", () => shouldBe(page, el("\"Edit chart\" dialog"), "visible"));
    await session.step(49, "And no point of the curve in row 1 of the \"Fitted Curve\" column should be left out", () => noPointLeftOut(page, 1, "Fitted Curve"));
    await session.step(50, "When user clicks on \"Close\" button in \"Edit chart\" dialog", () => clickOn(page, el("\"Close\" button in \"Edit chart\" dialog")));
    await session.step(51, "Then \"Edit chart\" dialog should be absent", () => shouldBe(page, el("\"Edit chart\" dialog"), "absent"));
    await session.step(52, "And no errors should have been logged", () => noErrors(page));
  });
});
