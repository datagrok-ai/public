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
import {bioInitialized, sequenceFontSize} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, shouldBe, uncheck, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {contextPanelClosed, openDatasetRows, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · To atomic level", () => {
  const session = feature(test, "features/demos/06-to-atomic-level.feature", import.meta.url);
  test("Sequences become molecules, monomer by monomer", {tag: ["@guide", "@demo", "@viewport:1440x810"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(18, "And user opens HELM_sample dataset keeping the first 5 rows", () => openDatasetRows(page, ds("HELM_sample"), 5));
    await session.step(19, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And the \"HELM\" column of the grid is 560 px wide", () => gridColumnWidth(page, "HELM", 560));
    await session.step(22, "When user picks \"Bio > Transform > To Atomic Level...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > To Atomic Level..."));
    await session.step(23, "And user unchecks \"Non-linear\" checkbox in \"To Atomic Level\" dialog", () => uncheck(page, el("\"Non-linear\" checkbox in \"To Atomic Level\" dialog")));
    await session.step(24, "And user checks \"Highlight monomers\" checkbox in \"To Atomic Level\" dialog", () => check(page, el("\"Highlight monomers\" checkbox in \"To Atomic Level\" dialog")));
    await session.step(25, "And user clicks on OK button in \"To Atomic Level\" dialog", () => clickOn(page, el("OK button in \"To Atomic Level\" dialog")));
    await session.step(26, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(27, "And a new column \"molfile(HELM)\" should have been added", () => newColumnNamed(page, "molfile(HELM)"));
    await session.step(28, "And \"molfile(HELM)\" column should have no missing values", () => columnComplete(page, "molfile(HELM)"));
    await session.step(29, "Given the grid shows only the columns \"molfile(HELM), HELM\"", () => gridShowsColumns(page, "molfile(HELM), HELM"));
    await session.step(30, "And the \"molfile(HELM)\" column of the grid is 720 px wide", () => gridColumnWidth(page, "molfile(HELM)", 720));
    await session.step(31, "And the \"HELM\" column of the grid is 560 px wide", () => gridColumnWidth(page, "HELM", 560));
    await session.step(32, "And the rows of the grid are 230 px tall", () => gridRowHeight(page, 230));
    await session.step(33, "And no row of the table is current", () => noCurrentRow(page));
    await session.step(34, "When the video says \"Every atom colored by the monomer it comes from\"", () => videoSays(page, "Every atom colored by the monomer it comes from"));
    await session.step(35, "And user hovers over monomer 1 of the \"cell 2 of HELM\" area of grid", () => hoverMonomer(page, 1, "cell 2 of HELM", el("grid")));
    await session.step(36, "And user hovers over monomer 7 of the \"cell 2 of HELM\" area of grid", () => hoverMonomer(page, 7, "cell 2 of HELM", el("grid")));
    await session.step(37, "And user hovers over monomer 14 of the \"cell 2 of HELM\" area of grid", () => hoverMonomer(page, 14, "cell 2 of HELM", el("grid")));
    await session.step(38, "And the video says \"Non-linear mode: cycles and bridges drawn as rings\"", () => videoSays(page, "Non-linear mode: cycles and bridges drawn as rings"));
    await session.step(39, "Given user opens BILN_complex dataset keeping the first 6 rows", () => openDatasetRows(page, ds("BILN_complex"), 6));
    await session.step(40, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(41, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(42, "And the grid shows only the columns \"name, biln\"", () => gridShowsColumns(page, "name, biln"));
    await session.step(43, "And the monomers of the \"biln\" column are drawn at 16 px", () => sequenceFontSize(page, "biln", 16));
    await session.step(44, "And the \"biln\" column of the grid is 560 px wide", () => gridColumnWidth(page, "biln", 560));
    await session.step(45, "And the rows of the grid are 60 px tall", () => gridRowHeight(page, 60));
    await session.step(46, "When user picks \"Bio > Transform > To Atomic Level...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > To Atomic Level..."));
    await session.step(47, "Then \"Non-linear\" checkbox in \"To Atomic Level\" dialog should be checked", () => shouldBe(page, el("\"Non-linear\" checkbox in \"To Atomic Level\" dialog"), "checked"));
    await session.step(48, "When user checks \"Highlight monomers\" checkbox in \"To Atomic Level\" dialog", () => check(page, el("\"Highlight monomers\" checkbox in \"To Atomic Level\" dialog")));
    await session.step(49, "And user clicks on OK button in \"To Atomic Level\" dialog", () => clickOn(page, el("OK button in \"To Atomic Level\" dialog")));
    await session.step(50, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(51, "And a new column \"molfile(biln)\" should have been added", () => newColumnNamed(page, "molfile(biln)"));
    await session.step(52, "And \"molfile(biln)\" column should have no missing values", () => columnComplete(page, "molfile(biln)"));
    await session.step(53, "Given the grid shows only the columns \"name, molfile(biln), biln\"", () => gridShowsColumns(page, "name, molfile(biln), biln"));
    await session.step(54, "And the \"molfile(biln)\" column of the grid is 560 px wide", () => gridColumnWidth(page, "molfile(biln)", 560));
    await session.step(55, "And the \"biln\" column of the grid is 560 px wide", () => gridColumnWidth(page, "biln", 560));
    await session.step(56, "And the rows of the grid are 230 px tall", () => gridRowHeight(page, 230));
    await session.step(57, "And no row of the table is current", () => noCurrentRow(page));
    await session.step(58, "When user hovers over monomer 1 of the \"cell 1 of biln\" area of grid", () => hoverMonomer(page, 1, "cell 1 of biln", el("grid")));
    await session.step(59, "And user hovers over monomer 4 of the \"cell 2 of biln\" area of grid", () => hoverMonomer(page, 4, "cell 2 of biln", el("grid")));
    await session.step(60, "Then no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
