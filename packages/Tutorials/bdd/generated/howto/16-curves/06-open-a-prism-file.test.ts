/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/16-curves/06-open-a-prism-file.feature
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
import {clickBesideCurve, doubleClickBesideCurve} from '../../../bindings/howto-curves.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, expand, isExpanded, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelOpen, contextPanelClosed, packageInstalled, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Curves", () => {
  const session = feature(test, "features/howto/16-curves/06-open-a-prism-file.feature", import.meta.url);
  test("Open a GraphPad Prism file", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(21, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(24, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(26, "When user expands Files---App-Data---Curves tree node inside browse tree", () => expand(page, el("Files---App-Data---Curves tree node inside browse tree")), undefined, "Prism files (.pzfx) open like any other file");
    await session.step(27, "And user double-clicks on Files---App-Data---Curves---dose-response-4pl.pzfx tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Curves---dose-response-4pl.pzfx tree node inside browse tree")));
    await session.step(28, "Then the \"dose-response-4pl.pzfx\" view should be current", () => viewIsCurrent(page, "dose-response-4pl.pzfx"));
    await session.step(29, "When the video says \"Each XY table is a row, its data sets fitted as curves\"", () => videoSays(page, "Each XY table is a row, its data sets fitted as curves"));
    await session.step(31, "And user double-clicks on the \"cell 1 of Fitted Curve\" area of grid beside its curve", () => doubleClickBesideCurve(page, "cell 1 of Fitted Curve", el("grid")), undefined, "Double-click the cell to see the curves larger");
    await session.step(32, "Then \"Edit chart\" dialog should be visible", () => shouldBe(page, el("\"Edit chart\" dialog"), "visible"));
    await session.step(33, "When the video says \"One curve per compound; values excluded in Prism are crosses\"", () => videoSays(page, "One curve per compound; values excluded in Prism are crosses"));
    await session.step(34, "And user clicks on \"Close\" button in \"Edit chart\" dialog", () => clickOn(page, el("\"Close\" button in \"Edit chart\" dialog")));
    await session.step(35, "Then \"Edit chart\" dialog should be absent", () => shouldBe(page, el("\"Edit chart\" dialog"), "absent"));
    await session.step(37, "When user clicks on the \"cell 1 of Fitted Curve\" area of grid beside its curve", () => clickBesideCurve(page, "cell 1 of Fitted Curve", el("grid")), undefined, "Click the cell, then open the context panel");
    await session.step(38, "And user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")));
    await session.step(39, "Then \"Fit\" pane in context panel should be visible", () => shouldBe(page, el("\"Fit\" pane in context panel"), "visible"));
    await session.step(41, "When user expands Fit accordion header in context panel", () => expand(page, el("Fit accordion header in context panel")), undefined, "Fit: R², IC50, slope and more for every compound");
    await session.step(42, "Then \"IC50\" input in context panel should be visible", () => shouldBe(page, el("\"IC50\" input in context panel"), "visible"));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
  });
});
