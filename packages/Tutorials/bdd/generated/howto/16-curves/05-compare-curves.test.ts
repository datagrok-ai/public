/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/16-curves/05-compare-curves.feature
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
import {gridScrolledToRow} from '../../../bindings/howto-curves.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, typeInto, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {onlyOfSelected, selectRowsByHeaders, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, contextPanelClosed, openDataset, packageInstalled, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noErrors, readingIs, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Curves", () => {
  const session = feature(test, "features/howto/16-curves/05-compare-curves.feature", import.meta.url);
  test("Compare curves in one chart", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(21, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(22, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(23, "And user opens assay-curves dataset", () => openDataset(page, ds("assay-curves")));
    await session.step(24, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(25, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(26, "And the grid shows only the columns \"Molecule, Fitted Curve, Compound ID, Run ID\"", () => gridShowsColumns(page, "Molecule, Fitted Curve, Compound ID, Run ID"));
    await session.step(27, "And the grid of the current view starts at row 34", () => gridScrolledToRow(page, 34));
    await session.step(29, "When user selects rows 34 to 37 by their row headers", () => selectRowsByHeaders(page, 34, 37), undefined, "Select the curves to compare: four runs of compound_id 20");
    await session.step(30, "Then 4 rows should be selected", () => selectedRowCount(page, 4));
    await session.step(31, "And only rows where \"Compound ID\" is \"compound_id 20\" should be selected", () => onlyOfSelected(page, "Compound ID", "compound_id 20"));
    await session.step(33, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")), undefined, "Add viewer, then search for curve");
    await session.step(34, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(35, "When user types \"curve\" into viewer gallery search", () => typeInto(page, "curve", el("viewer gallery search")));
    await session.step(37, "And user clicks on first \"MultiCurveViewer\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"MultiCurveViewer\" card in \"Add Viewer\" dialog")), undefined, "Multi curve viewer: many curves in one chart");
    await session.step(38, "Then \"Add Viewer\" dialog should be absent", () => shouldBe(page, el("\"Add Viewer\" dialog"), "absent"));
    await session.step(39, "And the open tableview should have 1 MultiCurveViewer viewer", () => viewerCount(page, 1, "MultiCurveViewer"));
    await session.step(40, "And the \"curves\" reading of MultiCurveViewer viewer should be 4", () => readingIs(page, "curves", el("MultiCurveViewer viewer"), 4));
    await session.step(41, "When the video says \"The four runs overlaid, the legend names each run\" at the bottom", () => videoSaysAt(page, "The four runs overlaid, the legend names each run", "bottom"));
    await session.step(42, "And the video says \"Run 3 (squares) sits left of the others: a lower IC50\" at the bottom", () => videoSaysAt(page, "Run 3 (squares) sits left of the others: a lower IC50", "bottom"));
    await session.step(44, "And user hovers over the \"cell 38 of Compound ID\" area of grid", () => hoverArea(page, "cell 38 of Compound ID", el("grid")), undefined, "Point at another row: its curve joins the comparison");
    await session.step(45, "Then the \"curves\" reading of MultiCurveViewer viewer should be 5", () => readingIs(page, "curves", el("MultiCurveViewer viewer"), 5));
    await session.step(46, "When the video says \"In red: compound_id 21, against compound_id 20's runs\" at the bottom", () => videoSaysAt(page, "In red: compound_id 21, against compound_id 20's runs", "bottom"));
    await session.step(47, "Then no errors should have been logged", () => noErrors(page));
  });
});
