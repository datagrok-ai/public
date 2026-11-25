/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/06-to-atomic-level.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '../../bindings/monomer-form.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {noCurrentRow} from '../../bindings/demo-atomic.js';
import {monomerNameLength} from '../../bindings/demo-notation.js';
import {bioInitialized, sequenceFontSize} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, shouldBe, uncheck, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, setGridColumnFormat} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelClosed, contextPanelClosed, openDatasetRows, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, pointerAway} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("2 · Every sequence is a molecule · To atomic level", () => {
  const session = feature(test, "features/demos/06-to-atomic-level.feature", import.meta.url);
  test("Sequences become molecules, monomer by monomer", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(19, "And user opens HELM_sample dataset keeping the first 5 rows", () => openDatasetRows(page, ds("HELM_sample"), 5));
    await session.step(20, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(23, "And the \"HELM\" column of the grid is 560 px wide", () => gridColumnWidth(page, "HELM", 560));
    await session.step(24, "And user sets the format of \"Activity\" column in the grid to \"0.00\"", () => setGridColumnFormat(page, "Activity", "0.00"));
    await session.step(25, "When user picks \"Bio > Transform > To Atomic Level...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > To Atomic Level..."));
    await session.step(26, "And user unchecks \"Non-linear\" checkbox in \"To Atomic Level\" dialog", () => uncheck(page, el("\"Non-linear\" checkbox in \"To Atomic Level\" dialog")));
    await session.step(27, "And user checks \"Highlight monomers\" checkbox in \"To Atomic Level\" dialog", () => check(page, el("\"Highlight monomers\" checkbox in \"To Atomic Level\" dialog")));
    await session.step(28, "And user clicks on OK button in \"To Atomic Level\" dialog", () => clickOn(page, el("OK button in \"To Atomic Level\" dialog")));
    await session.step(29, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(30, "And a new column \"molfile(HELM)\" should have been added", () => newColumnNamed(page, "molfile(HELM)"));
    await session.step(31, "And \"molfile(HELM)\" column should have no missing values", () => columnComplete(page, "molfile(HELM)"));
    await session.step(32, "Given the grid shows only the columns \"molfile(HELM), HELM\"", () => gridShowsColumns(page, "molfile(HELM), HELM"));
    await session.step(33, "And the \"molfile(HELM)\" column of the grid is 820 px wide", () => gridColumnWidth(page, "molfile(HELM)", 820));
    await session.step(34, "And the \"HELM\" column of the grid is 620 px wide", () => gridColumnWidth(page, "HELM", 620));
    await session.step(35, "And the rows of the grid are 230 px tall", () => gridRowHeight(page, 230));
    await session.step(36, "And no row of the table is current", () => noCurrentRow(page));
    await session.step(37, "When the video says \"Every atom colored by the monomer it comes from\"", () => videoSays(page, "Every atom colored by the monomer it comes from"));
    await session.step(38, "And user hovers over monomer 1 of the \"cell 2 of HELM\" area of grid", () => hoverMonomer(page, 1, "cell 2 of HELM", el("grid")));
    await session.step(39, "And user hovers over monomer 7 of the \"cell 2 of HELM\" area of grid", () => hoverMonomer(page, 7, "cell 2 of HELM", el("grid")));
    await session.step(40, "And user hovers over monomer 14 of the \"cell 2 of HELM\" area of grid", () => hoverMonomer(page, 14, "cell 2 of HELM", el("grid")));
    await session.step(41, "Given user opens BILN_complex dataset keeping the first 6 rows", () => openDatasetRows(page, ds("BILN_complex"), 6));
    await session.step(42, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(43, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(44, "And the grid shows only the columns \"name, biln\"", () => gridShowsColumns(page, "name, biln"));
    await session.step(45, "And the monomers of the \"biln\" column are drawn at 16 px", () => sequenceFontSize(page, "biln", 16));
    await session.step(46, "And the monomer names of the \"biln\" column are cut at 8 characters", () => monomerNameLength(page, "biln", 8));
    await session.step(47, "And the \"biln\" column of the grid is 560 px wide", () => gridColumnWidth(page, "biln", 560));
    await session.step(48, "And the rows of the grid are 60 px tall", () => gridRowHeight(page, 60));
    await session.step(49, "And user picks \"Bio > Transform > To Atomic Level...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > To Atomic Level..."));
    await session.step(50, "Then \"Non-linear\" checkbox in \"To Atomic Level\" dialog should be checked", () => shouldBe(page, el("\"Non-linear\" checkbox in \"To Atomic Level\" dialog"), "checked"));
    await session.step(51, "Given user checks \"Highlight monomers\" checkbox in \"To Atomic Level\" dialog", () => check(page, el("\"Highlight monomers\" checkbox in \"To Atomic Level\" dialog")));
    await session.step(52, "And user clicks on OK button in \"To Atomic Level\" dialog", () => clickOn(page, el("OK button in \"To Atomic Level\" dialog")));
    await session.step(53, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(54, "And a new column \"molfile(biln)\" should have been added", () => newColumnNamed(page, "molfile(biln)"));
    await session.step(55, "And \"molfile(biln)\" column should have no missing values", () => columnComplete(page, "molfile(biln)"));
    await session.step(56, "Given the grid shows only the columns \"name, molfile(biln), biln\"", () => gridShowsColumns(page, "name, molfile(biln), biln"));
    await session.step(57, "And the \"molfile(biln)\" column of the grid is 640 px wide", () => gridColumnWidth(page, "molfile(biln)", 640));
    await session.step(58, "And the \"biln\" column of the grid is 600 px wide", () => gridColumnWidth(page, "biln", 600));
    await session.step(59, "And the rows of the grid are 230 px tall", () => gridRowHeight(page, 230));
    await session.step(60, "And no row of the table is current", () => noCurrentRow(page));
    await session.step(61, "When the video says \"Non-linear mode: cycles and bridges drawn as rings\"", () => videoSays(page, "Non-linear mode: cycles and bridges drawn as rings"));
    await session.step(62, "And user hovers over monomer 1 of the \"cell 1 of biln\" area of grid", () => hoverMonomer(page, 1, "cell 1 of biln", el("grid")));
    await session.step(63, "And user hovers over monomer 4 of the \"cell 2 of biln\" area of grid", () => hoverMonomer(page, 4, "cell 2 of biln", el("grid")));
    await session.step(64, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(65, "Then no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
