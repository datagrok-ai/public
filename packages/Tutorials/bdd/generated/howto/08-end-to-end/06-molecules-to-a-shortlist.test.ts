/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/08-end-to-end/06-molecules-to-a-shortlist.feature
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
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, pressKeyIn, shouldBe, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPasses, rowCount, selectedPassFilter, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, openDataset, packageInstalled, simpleModeOff, sketcherIs, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dragSelectionBetweenAreas, hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · End to end", () => {
  const session = feature(test, "features/howto/08-end-to-end/06-molecules-to-a-shortlist.feature", import.meta.url);
  test("From molecules to a shortlist", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(19, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(20, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(21, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(22, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(23, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(24, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(25, "When the video says \"The task: a shortlist of potent compounds with a pyrrolidine ring\"", () => videoSays(page, "The task: a shortlist of potent compounds with a pyrrolidine ring"));
    await session.step(27, "And user picks \"Chem > Search > Substructure Search...\" from the top menu", () => pickFromTopMenu(page, "Chem > Search > Substructure Search..."), undefined, "Keep the molecules with the ring: Chem > Search > Substructure Search");
    await session.step(28, "Then \"Substructure search\" dialog should be visible", () => shouldBe(page, el("\"Substructure search\" dialog"), "visible"));
    await session.step(29, "When user clicks on OK button in \"Substructure search\" dialog", () => clickOn(page, el("OK button in \"Substructure search\" dialog")));
    await session.step(30, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(31, "When user clicks on molecule input of sketcher dialog", () => clickOn(page, el("molecule input of sketcher dialog")));
    await session.step(32, "And user types \"C1CCNC1\" into molecule input of sketcher dialog", () => typeInto(page, "C1CCNC1", el("molecule input of sketcher dialog")));
    await session.step(33, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(34, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(35, "Then 33 rows should pass the filter", () => filterPasses(page, 33));
    await session.step(36, "When the video says \"33 of 100 have the ring\"", () => videoSays(page, "33 of 100 have the ring"));
    await session.step(38, "And user clicks on \"scatter plot\" icon in toolbox", () => clickOn(page, el("\"scatter plot\" icon in toolbox")), undefined, "Plot two properties: logP against polar surface area");
    await session.step(39, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(40, "When user picks \"NIBR logP\" in the \"x\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "NIBR logP", "x", el("scatter plot viewer")));
    await session.step(41, "And user picks \"TPSA\" in the \"y\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "TPSA", "y", el("scatter plot viewer")));
    await session.step(43, "And user picks \"Competition assay\" in the \"color\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "Competition assay", "color", el("scatter plot viewer")), undefined, "Color by activity. Lower Competition assay values are more potent");
    await session.step(45, "And user clicks on \"histogram\" icon in toolbox", () => clickOn(page, el("\"histogram\" icon in toolbox")), undefined, "A histogram of the activity: the most potent are on the left");
    await session.step(46, "Then histogram viewer should be visible", () => shouldBe(page, el("histogram viewer"), "visible"));
    await session.step(47, "When user picks \"Competition assay\" in the \"value\" column selector of histogram viewer", () => pickInColumnSelector(page, "Competition assay", "value", el("histogram viewer")));
    await session.step(49, "And user drags a selection box from the \"bin 1\" area to the \"bin 2\" area of histogram viewer", () => dragSelectionBetweenAreas(page, "bin 1", "bin 2", el("histogram viewer")), undefined, "Shift+drag over the first bars to select the most potent");
    await session.step(50, "Then 7 rows should be selected", () => selectedRowCount(page, 7));
    await session.step(51, "And every selected row should pass the filter", () => selectedPassFilter(page));
    await session.step(52, "When the video says \"The plot marks where they sit\"", () => videoSays(page, "The plot marks where they sit"));
    await session.step(53, "And user hovers over scatter plot viewer", () => hoverOver(page, el("scatter plot viewer")));
    await session.step(55, "And user picks \"Select > Extract Selected Rows\" from the top menu", () => pickFromTopMenu(page, "Select > Extract Selected Rows"), undefined, "Select > Extract Selected Rows: the shortlist as a table of its own");
    await session.step(56, "Then the table should have 7 rows", () => rowCount(page, 7));
    await session.step(57, "When the video says \"The shortlist, with structures and every property\"", () => videoSays(page, "The shortlist, with structures and every property"));
    await session.step(58, "And user hovers over the \"cell 1 of Structure\" area of grid", () => hoverArea(page, "cell 1 of Structure", el("grid")));
    await session.step(59, "Then no errors should have been logged", () => noErrors(page));
  });
});
