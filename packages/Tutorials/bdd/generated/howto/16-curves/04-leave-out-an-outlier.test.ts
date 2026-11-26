/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/16-curves/04-leave-out-an-outlier.feature
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
import {clickCurvePoint, gridScrolledToRow, hoverCurvePoint, pointKept, pointLeftOut, statisticColumnAdded} from '../../../bindings/howto-curves.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {displayedInRow, setGridColumnFormat} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, contextPanelClosed, openDataset, packageInstalled, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Curves", () => {
  const session = feature(test, "features/howto/16-curves/04-leave-out-an-outlier.feature", import.meta.url);
  test("Leave out an outlier", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(20, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(21, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(22, "And user opens assay-curves dataset", () => openDataset(page, ds("assay-curves")));
    await session.step(23, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(24, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(25, "And the \"ic50\" statistic of the \"Fitted Curve\" column is a column of the table", () => statisticColumnAdded(page, "ic50", "Fitted Curve"));
    await session.step(26, "And the \"rSquared\" statistic of the \"Fitted Curve\" column is a column of the table", () => statisticColumnAdded(page, "rSquared", "Fitted Curve"));
    await session.step(27, "And user sets the format of \"Fitted Curve 1 ic50\" column in the grid to \"scientific\"", () => setGridColumnFormat(page, "Fitted Curve 1 ic50", "scientific"));
    await session.step(28, "And the grid shows only the columns \"Molecule, Fitted Curve, Fitted Curve 1 rSquared, Fitted Curve 1 ic50, Compound|Assay|Target\"", () => gridShowsColumns(page, "Molecule, Fitted Curve, Fitted Curve 1 rSquared, Fitted Curve 1 ic50, Compound|Assay|Target"));
    await session.step(29, "And the rows of the grid are 230 px tall", () => gridRowHeight(page, 230));
    await session.step(30, "And the \"Fitted Curve\" column of the grid is 340 px wide", () => gridColumnWidth(page, "Fitted Curve", 340));
    await session.step(31, "And the grid of the current view starts at row 24", () => gridScrolledToRow(page, 24));
    await session.step(32, "Then the \"Fitted Curve 1 rSquared\" cell of row 24 should be displayed as \"0.93\"", () => displayedInRow(page, "Fitted Curve 1 rSquared", 24, "0.93"));
    await session.step(33, "And the \"Fitted Curve 1 ic50\" cell of row 24 should be displayed as \"4.96E-8\"", () => displayedInRow(page, "Fitted Curve 1 ic50", 24, "4.96E-8"));
    await session.step(34, "When the video says \"R² and IC50 of every curve, calculated from its fit\"", () => videoSays(page, "R² and IC50 of every curve, calculated from its fit"));
    await session.step(36, "And user hovers over the \"cell 24 of Fitted Curve 1 rSquared\" area of grid", () => hoverArea(page, "cell 24 of Fitted Curve 1 rSquared", el("grid")), undefined, "R² 0.93: this curve misses one of its points");
    await session.step(38, "And user hovers over point 5 of the curve in the \"cell 24 of Fitted Curve\" area of grid", () => hoverCurvePoint(page, 5, "cell 24 of Fitted Curve", el("grid")), undefined, "This point sits well below its neighbours");
    await session.step(40, "And user clicks on point 5 of the curve in the \"cell 24 of Fitted Curve\" area of grid", () => clickCurvePoint(page, 5, "cell 24 of Fitted Curve", el("grid")), undefined, "Click it to leave it out of the fit");
    await session.step(41, "Then point 5 of the curve in row 24 of the \"Fitted Curve\" column should be left out", () => pointLeftOut(page, 5, 24, "Fitted Curve"));
    await session.step(42, "And the \"Fitted Curve 1 rSquared\" cell of row 24 should be displayed as \"0.98\"", () => displayedInRow(page, "Fitted Curve 1 rSquared", 24, "0.98"));
    await session.step(43, "And the \"Fitted Curve 1 ic50\" cell of row 24 should be displayed as \"4.25E-8\"", () => displayedInRow(page, "Fitted Curve 1 ic50", 24, "4.25E-8"));
    await session.step(44, "When the video says \"Left out: a red cross, and the curve is fitted without it\"", () => videoSays(page, "Left out: a red cross, and the curve is fitted without it"));
    await session.step(46, "And user hovers over the \"cell 24 of Fitted Curve 1 rSquared\" area of grid", () => hoverArea(page, "cell 24 of Fitted Curve 1 rSquared", el("grid")), undefined, "R² follows: 0.98 now");
    await session.step(48, "And user hovers over the \"cell 24 of Fitted Curve 1 ic50\" area of grid", () => hoverArea(page, "cell 24 of Fitted Curve 1 ic50", el("grid")), undefined, "So does IC50: from 4.96E-8 to 4.25E-8");
    await session.step(50, "And user clicks on point 5 of the curve in the \"cell 24 of Fitted Curve\" area of grid", () => clickCurvePoint(page, 5, "cell 24 of Fitted Curve", el("grid")), undefined, "Click the cross to put the point back");
    await session.step(51, "Then point 5 of the curve in row 24 of the \"Fitted Curve\" column should be in the fit", () => pointKept(page, 5, 24, "Fitted Curve"));
    await session.step(52, "And the \"Fitted Curve 1 rSquared\" cell of row 24 should be displayed as \"0.93\"", () => displayedInRow(page, "Fitted Curve 1 rSquared", 24, "0.93"));
    await session.step(53, "And the \"Fitted Curve 1 ic50\" cell of row 24 should be displayed as \"4.96E-8\"", () => displayedInRow(page, "Fitted Curve 1 ic50", 24, "4.96E-8"));
    await session.step(54, "When the video says \"Back in the fit: R² 0.93 and IC50 4.96E-8 again\"", () => videoSays(page, "Back in the fit: R² 0.93 and IC50 4.96E-8 again"));
    await session.step(55, "Then no errors should have been logged", () => noErrors(page));
  });
});
