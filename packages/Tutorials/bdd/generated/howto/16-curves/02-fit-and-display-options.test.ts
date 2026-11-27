/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/16-curves/02-fit-and-display-options.feature
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
import {choiceTicked, clickBesideCurve, noPointLeftOut, tickChoice} from '../../../bindings/howto-curves.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, expand, selectIn, shouldBe, shouldHaveValue, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, contextPanelOpen, openDataset, packageInstalled, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Curves", () => {
  const session = feature(test, "features/howto/16-curves/02-fit-and-display-options.feature", import.meta.url);
  test("Change how curves are fitted and drawn", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(18, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And user opens assay-curves dataset", () => openDataset(page, ds("assay-curves")));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the rows of the grid are 230 px tall", () => gridRowHeight(page, 230));
    await session.step(23, "And the \"Fitted Curve\" column of the grid is 340 px wide", () => gridColumnWidth(page, "Fitted Curve", 340));
    await session.step(24, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(26, "When user clicks on the \"cell 1 of Fitted Curve\" area of grid", () => clickArea(page, "cell 1 of Fitted Curve", el("grid")), undefined, "Click a curve: its options are in the context panel");
    await session.step(27, "Then no point of the curve in row 1 of the \"Fitted Curve\" column should be left out", () => noPointLeftOut(page, 1, "Fitted Curve"));
    await session.step(28, "And \"Options\" pane in context panel should be visible", () => shouldBe(page, el("\"Options\" pane in context panel"), "visible"));
    await session.step(29, "When user expands Options accordion header in context panel", () => expand(page, el("Options accordion header in context panel")));
    await session.step(30, "Then \"Level\" input in context panel should have value \"Column\"", () => shouldHaveValue(page, el("\"Level\" input in context panel"), "Column"));
    await session.step(31, "When the video says \"Level: Column, so each change reaches every curve in the column\"", () => videoSays(page, "Level: Column, so each change reaches every curve in the column"));
    await session.step(33, "And user selects \"linear\" in \"Fit Function\" input in context panel", () => selectIn(page, "linear", el("\"Fit Function\" input in context panel")), undefined, "Fit Function: the model fitted to the points");
    await session.step(34, "Then \"Fit Function\" input in context panel should have value \"linear\"", () => shouldHaveValue(page, el("\"Fit Function\" input in context panel"), "linear"));
    await session.step(36, "When user selects \"4pl-dose-response\" in \"Fit Function\" input in context panel", () => selectIn(page, "4pl-dose-response", el("\"Fit Function\" input in context panel")), undefined, "Back to the four-parameter dose-response fit");
    await session.step(38, "And user checks \"Confidence Interval\" input in context panel", () => check(page, el("\"Confidence Interval\" input in context panel")), undefined, "Confidence Interval: a band around each fit");
    await session.step(40, "And user ticks \"IC80\" in \"Droplines\" input in context panel", () => tickChoice(page, "IC80", el("\"Droplines\" input in context panel")), undefined, "Droplines: IC80 next to IC50");
    await session.step(42, "And user ticks \"ic50\" in \"Statistics\" input in context panel", () => tickChoice(page, "ic50", el("\"Statistics\" input in context panel")), undefined, "Statistics: values printed on every plot");
    await session.step(43, "And user ticks \"rSquared\" in \"Statistics\" input in context panel", () => tickChoice(page, "rSquared", el("\"Statistics\" input in context panel")));
    await session.step(45, "And user clicks on the \"cell 2 of Fitted Curve\" area of grid beside its curve", () => clickBesideCurve(page, "cell 2 of Fitted Curve", el("grid")), undefined, "Another curve of the column has the same options");
    await session.step(46, "Then no point of the curve in row 2 of the \"Fitted Curve\" column should be left out", () => noPointLeftOut(page, 2, "Fitted Curve"));
    await session.step(47, "And \"Fit Function\" input in context panel should have value \"4pl-dose-response\"", () => shouldHaveValue(page, el("\"Fit Function\" input in context panel"), "4pl-dose-response"));
    await session.step(48, "And \"Confidence Interval\" input in context panel should be checked", () => shouldBe(page, el("\"Confidence Interval\" input in context panel"), "checked"));
    await session.step(49, "And \"ic50\" should be ticked in \"Statistics\" input in context panel", () => choiceTicked(page, "ic50", el("\"Statistics\" input in context panel")));
    await session.step(50, "And no errors should have been logged", () => noErrors(page));
  });
});
