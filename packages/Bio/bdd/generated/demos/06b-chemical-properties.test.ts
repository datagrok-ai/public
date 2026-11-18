/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/06b-chemical-properties.feature
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
import {atomicLevelDefaults, clickMarker} from '../../bindings/demo-chem.js';
import {monomerNameLength} from '../../bindings/demo-notation.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, shouldBe, shouldContainText, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, everyValueBetween, setGridColumnFormat} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {currentRowValue} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDatasetRows, packageInstalled, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noBalloons, pointerAway, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("2 · Every sequence is a molecule · Chemical properties", () => {
  const session = feature(test, "features/demos/06b-chemical-properties.feature", import.meta.url);
  test("Peptide drugs measured like small molecules", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(21, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(22, "And user opens BILN_complex dataset keeping the first 63 rows", () => openDatasetRows(page, ds("BILN_complex"), 63));
    await session.step(23, "And the \"biln\" column is converted to atomic level", () => atomicLevelDefaults(page, "biln"));
    await session.step(24, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(25, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(26, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(27, "And the grid shows only the columns \"name, category, biln, molfile(biln)\"", () => gridShowsColumns(page, "name, category, biln, molfile(biln)"));
    await session.step(28, "And the \"biln\" column of the grid is 520 px wide", () => gridColumnWidth(page, "biln", 520));
    await session.step(29, "And the \"molfile(biln)\" column of the grid is 360 px wide", () => gridColumnWidth(page, "molfile(biln)", 360));
    await session.step(30, "And the rows of the grid are 90 px tall", () => gridRowHeight(page, 90));
    await session.step(31, "When user picks \"Chem > Calculate > Chemical Properties...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Chemical Properties..."));
    await session.step(32, "Then \"Chemical Properties\" dialog should be visible", () => shouldBe(page, el("\"Chemical Properties\" dialog"), "visible"));
    await session.step(33, "And Molecules input in \"Chemical Properties\" dialog should contain text \"molfile(biln)\"", () => shouldContainText(page, el("Molecules input in \"Chemical Properties\" dialog"), "molfile(biln)"));
    await session.step(34, "When user checks \"Chemical Properties\" calculator in \"Chemical Properties\" dialog", () => check(page, el("\"Chemical Properties\" calculator in \"Chemical Properties\" dialog")));
    await session.step(35, "Then MW input in \"Chemical Properties\" dialog should be checked", () => shouldBe(page, el("MW input in \"Chemical Properties\" dialog"), "checked"));
    await session.step(36, "When user checks PSA input in \"Chemical Properties\" dialog", () => check(page, el("PSA input in \"Chemical Properties\" dialog")));
    await session.step(37, "And user checks \"Log P\" input in \"Chemical Properties\" dialog", () => check(page, el("\"Log P\" input in \"Chemical Properties\" dialog")));
    await session.step(38, "And user clicks on OK button in \"Chemical Properties\" dialog", () => clickOn(page, el("OK button in \"Chemical Properties\" dialog")));
    await session.step(39, "Then \"Chemical Properties\" dialog should be hidden", () => shouldBe(page, el("\"Chemical Properties\" dialog"), "hidden"));
    await session.step(40, "And a new column \"MW\" should have been added", () => newColumnNamed(page, "MW"));
    await session.step(41, "And a new column \"PSA\" should have been added", () => newColumnNamed(page, "PSA"));
    await session.step(42, "And a new column \"LogP\" should have been added", () => newColumnNamed(page, "LogP"));
    await session.step(43, "And \"MW\" column should have no missing values", () => columnComplete(page, "MW"));
    await session.step(44, "And every value of \"MW\" column should lie between 500 and 10000", () => everyValueBetween(page, "MW", 500, 10000));
    await session.step(45, "Given user sets the format of \"MW\" column in the grid to \"0\"", () => setGridColumnFormat(page, "MW", "0"));
    await session.step(46, "And user sets the format of \"PSA\" column in the grid to \"0\"", () => setGridColumnFormat(page, "PSA", "0"));
    await session.step(47, "And user sets the format of \"LogP\" column in the grid to \"0.0\"", () => setGridColumnFormat(page, "LogP", "0.0"));
    await session.step(48, "And the grid shows only the columns \"name, category, MW, PSA, LogP, biln\"", () => gridShowsColumns(page, "name, category, MW, PSA, LogP, biln"));
    await session.step(49, "And the monomer names of the \"biln\" column are cut at 8 characters", () => monomerNameLength(page, "biln", 8));
    await session.step(50, "When the video says \"Weight, polar surface and LogP for every peptide\" at the top", () => videoSaysAt(page, "Weight, polar surface and LogP for every peptide", "top"));
    await session.step(51, "Given the \"biln\" column of the grid is 300 px wide", () => gridColumnWidth(page, "biln", 300));
    await session.step(52, "And the rows of the grid are 40 px tall", () => gridRowHeight(page, 40));
    await session.step(53, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","MW"],["yColumnName","PSA"],["colorColumnName","category"],["xAxisType","logarithmic"],["yAxisType","logarithmic"],["xMin","300"],["yMin","100"],["legendPosition","Top"]]), [["xColumnName","MW"],["yColumnName","PSA"],["colorColumnName","category"],["xAxisType","logarithmic"],["yAxisType","logarithmic"],["xMin","300"],["yMin","100"],["legendPosition","Top"]]);
    await session.step(62, "And user sets \"formulaLines\" property of scatter plot viewer to '[{\"type\":\"line\",\"formula\":\"${MW} = 500\",\"title\":\"MW 500\",\"color\":\"#e15759\"},{\"type\":\"line\",\"formula\":\"${PSA} = 140\",\"title\":\"TPSA 140\",\"color\":\"#e15759\"}]'", () => setProperty(page, "formulaLines", el("scatter plot viewer"), "[{\"type\":\"line\",\"formula\":\"${MW} = 500\",\"title\":\"MW 500\",\"color\":\"#e15759\"},{\"type\":\"line\",\"formula\":\"${PSA} = 140\",\"title\":\"TPSA 140\",\"color\":\"#e15759\"}]"));
    await session.step(63, "And no row of the table is current", () => noCurrentRow(page));
    await session.step(64, "When user clicks on the marker of row 35 of scatter plot viewer", () => clickMarker(page, 35, el("scatter plot viewer")));
    await session.step(65, "Then \"name\" of the current row should be \"cyclosporin_like\"", () => currentRowValue(page, "name", "cyclosporin_like"));
    await session.step(66, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
    await session.step(67, "And the video says \"Cyclosporin and its peers: beyond the rule of five\" at the top", () => videoSaysAt(page, "Cyclosporin and its peers: beyond the rule of five", "top"));
    await session.step(68, "Then no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
