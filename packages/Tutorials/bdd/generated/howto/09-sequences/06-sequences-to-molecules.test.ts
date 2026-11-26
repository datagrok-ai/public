/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/09-sequences/06-sequences-to-molecules.feature
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
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {bioReady, noRowCurrent} from '../../../bindings/howto-sequences.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, shouldBe, shouldContainText, uncheck, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, setGridColumnFormat} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelClosed, contextPanelClosed, openDatasetRows, packageInstalled, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dragAreaBy, hoverArea, noBalloons, noErrors, pointerAway} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Sequences", () => {
  const session = feature(test, "features/howto/09-sequences/06-sequences-to-molecules.feature", import.meta.url);
  test("Turn sequences into molecules", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(17, "And the Bio package is initialized", () => bioReady(page));
    await session.step(18, "And user opens HELM_sample dataset keeping the first 5 rows", () => openDatasetRows(page, ds("HELM_sample"), 5));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(21, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(22, "And the \"HELM\" column of the grid is 520 px wide", () => gridColumnWidth(page, "HELM", 520));
    await session.step(23, "And user sets the format of \"Activity\" column in the grid to \"0.00\"", () => setGridColumnFormat(page, "Activity", "0.00"));
    await session.step(24, "And the rows of the grid are 150 px tall", () => gridRowHeight(page, 150));
    await session.step(25, "And no row of the table is current", () => noRowCurrent(page));
    await session.step(27, "When user picks \"Bio > Transform > To Atomic Level...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > To Atomic Level..."), undefined, "Bio > Transform > To Atomic Level");
    await session.step(28, "Then \"To Atomic Level\" dialog should be visible", () => shouldBe(page, el("\"To Atomic Level\" dialog"), "visible"));
    await session.step(30, "When user unchecks \"Non-linear\" checkbox in \"To Atomic Level\" dialog", () => uncheck(page, el("\"Non-linear\" checkbox in \"To Atomic Level\" dialog")), undefined, "Highlight monomers: color each atom by its monomer");
    await session.step(31, "And user checks \"Highlight monomers\" checkbox in \"To Atomic Level\" dialog", () => check(page, el("\"Highlight monomers\" checkbox in \"To Atomic Level\" dialog")));
    await session.step(32, "And user clicks on OK button in \"To Atomic Level\" dialog", () => clickOn(page, el("OK button in \"To Atomic Level\" dialog")));
    await session.step(33, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(34, "And a new column \"molfile(HELM)\" should have been added", () => newColumnNamed(page, "molfile(HELM)"));
    await session.step(35, "And \"molfile(HELM)\" column should have no missing values", () => columnComplete(page, "molfile(HELM)"));
    await session.step(37, "When user drags the \"column resizer molfile(HELM)\" area of grid by 420 pixels to the right", () => dragAreaBy(page, "column resizer molfile(HELM)", el("grid"), 420, "right"), undefined, "Drag the column's edge to see the molecules larger");
    await session.step(38, "And the video says \"Every atom colored by the monomer it comes from\"", () => videoSays(page, "Every atom colored by the monomer it comes from"));
    await session.step(39, "And user hovers over monomer 1 of the \"cell 2 of HELM\" area of grid", () => hoverMonomer(page, 1, "cell 2 of HELM", el("grid")));
    await session.step(40, "And user hovers over monomer 7 of the \"cell 2 of HELM\" area of grid", () => hoverMonomer(page, 7, "cell 2 of HELM", el("grid")));
    await session.step(41, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(42, "And the video says \"A molecule column now: the chemistry tools work on it\"", () => videoSays(page, "A molecule column now: the chemistry tools work on it"));
    await session.step(44, "And user picks \"Chem > Calculate > Chemical Properties...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Chemical Properties..."), undefined, "Chem > Calculate > Chemical Properties");
    await session.step(45, "Then \"Chemical Properties\" dialog should be visible", () => shouldBe(page, el("\"Chemical Properties\" dialog"), "visible"));
    await session.step(46, "And Molecules input in \"Chemical Properties\" dialog should contain text \"molfile(HELM)\"", () => shouldContainText(page, el("Molecules input in \"Chemical Properties\" dialog"), "molfile(HELM)"));
    await session.step(48, "When user checks \"Chemical Properties\" calculator in \"Chemical Properties\" dialog", () => check(page, el("\"Chemical Properties\" calculator in \"Chemical Properties\" dialog")), undefined, "Pick the properties: weight and LogP");
    await session.step(49, "Then MW input in \"Chemical Properties\" dialog should be checked", () => shouldBe(page, el("MW input in \"Chemical Properties\" dialog"), "checked"));
    await session.step(50, "When user checks \"Log P\" input in \"Chemical Properties\" dialog", () => check(page, el("\"Log P\" input in \"Chemical Properties\" dialog")));
    await session.step(51, "And user clicks on OK button in \"Chemical Properties\" dialog", () => clickOn(page, el("OK button in \"Chemical Properties\" dialog")));
    await session.step(52, "Then \"Chemical Properties\" dialog should be hidden", () => shouldBe(page, el("\"Chemical Properties\" dialog"), "hidden"));
    await session.step(53, "And a new column \"MW\" should have been added", () => newColumnNamed(page, "MW"));
    await session.step(54, "And a new column \"LogP\" should have been added", () => newColumnNamed(page, "LogP"));
    await session.step(55, "And \"MW\" column should have no missing values", () => columnComplete(page, "MW"));
    await session.step(56, "When the video says \"Weight and LogP for every peptide, like any small molecule\"", () => videoSays(page, "Weight and LogP for every peptide, like any small molecule"));
    await session.step(57, "And user hovers over the \"cell 2 of MW\" area of grid", () => hoverArea(page, "cell 2 of MW", el("grid")));
    await session.step(58, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(59, "And no errors should have been logged", () => noErrors(page));
  });
});
