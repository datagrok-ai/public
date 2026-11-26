/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/09-sequences/10-scan-liabilities.feature
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
import {bioReady, columnStaysWide} from '../../../bindings/howto-sequences.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {balloonsClosed, contextPanelClosed, openDataset, packageInstalled, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dragAreaBy, noBalloons, noErrors, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridRowHeight, gridShowsColumns, headerStartsAt, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Sequences", () => {
  const session = feature(test, "features/howto/09-sequences/10-scan-liabilities.feature", import.meta.url);
  test("Scan sequences for liabilities", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the \"HMM\" package is installed", () => packageInstalled(page, "HMM"));
    await session.step(17, "And user opens antibodies dataset", () => openDataset(page, ds("antibodies")));
    await session.step(18, "And the Bio package is initialized", () => bioReady(page));
    await session.step(19, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
    await session.step(22, "And user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
    await session.step(23, "Then a new column \"AntibodyHC (aligned)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)"));
    await session.step(24, "Given the notification balloons are closed", () => balloonsClosed(page));
    await session.step(25, "And the grid shows only the columns \"Antibody_ID, AntibodyHC (aligned)\"", () => gridShowsColumns(page, "Antibody_ID, AntibodyHC (aligned)"));
    await session.step(26, "And the header of the \"AntibodyHC (aligned)\" column starts at position \"95\"", () => headerStartsAt(page, "AntibodyHC (aligned)", "95"));
    await session.step(27, "And the \"AntibodyHC (aligned)\" column of the grid stays 560 px wide", () => columnStaysWide(page, "AntibodyHC (aligned)", 560));
    await session.step(28, "And the rows of the grid are 30 px tall", () => gridRowHeight(page, 30));
    await session.step(29, "When the video says \"Antibody heavy chains, numbered in IMGT\"", () => videoSays(page, "Antibody heavy chains, numbered in IMGT"));
    await session.step(31, "And user picks \"Annotations > Extract CDR3 as Column\" from the context menu of the \"cell 12 of AntibodyHC (aligned)\" area of grid", () => pickFromAreaContextMenu(page, "Annotations > Extract CDR3 as Column", "cell 12 of AntibodyHC (aligned)", el("grid")), undefined, "Right-click a chain: extract CDR3 as a column");
    await session.step(32, "Then a new column \"AntibodyHC (aligned)(CDR3)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)(CDR3)"));
    await session.step(34, "When user drags the \"column resizer AntibodyHC (aligned)(CDR3)\" area of grid by 300 pixels to the right", () => dragAreaBy(page, "column resizer AntibodyHC (aligned)(CDR3)", el("grid"), 300, "right"), undefined, "Drag its edge to widen it");
    await session.step(36, "When user picks \"Bio > Annotate > Scan Liabilities...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Scan Liabilities..."), undefined, "Bio > Annotate > Scan Liabilities");
    await session.step(37, "Then \"Scan Sequence Liabilities\" dialog should be visible", () => shouldBe(page, el("\"Scan Sequence Liabilities\" dialog"), "visible"));
    await session.step(39, "When user selects \"AntibodyHC (aligned)(CDR3)\" in Sequence input in \"Scan Sequence Liabilities\" dialog", () => selectIn(page, "AntibodyHC (aligned)(CDR3)", el("Sequence input in \"Scan Sequence Liabilities\" dialog")), undefined, "Pick the column to scan: the CDR3s");
    await session.step(40, "And user clicks on OK button in \"Scan Sequence Liabilities\" dialog", () => clickOn(page, el("OK button in \"Scan Sequence Liabilities\" dialog")));
    await session.step(41, "Then a new column \"~AntibodyHC (aligned)(CDR3)_annotations\" should have been added", () => newColumnNamed(page, "~AntibodyHC (aligned)(CDR3)_annotations"));
    await session.step(42, "When the video says \"Each liability found is underlined on its CDR3\"", () => videoSays(page, "Each liability found is underlined on its CDR3"));
    await session.step(43, "And user hovers over monomer 3 of the \"cell 12 of AntibodyHC (aligned)(CDR3)\" area of grid", () => hoverMonomer(page, 3, "cell 12 of AntibodyHC (aligned)(CDR3)", el("grid")));
    await session.step(44, "Then tooltip should contain text \"Isomerization (DG)\"", () => shouldContainText(page, el("tooltip"), "Isomerization (DG)"));
    await session.step(45, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
  });
});
