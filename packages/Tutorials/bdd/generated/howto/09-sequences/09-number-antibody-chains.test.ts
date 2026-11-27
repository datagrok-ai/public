/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/09-sequences/09-number-antibody-chains.feature
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
import {bioReady, hoverRegionBand} from '../../../bindings/howto-sequences.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveValue, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {contextPanelClosed, openDataset, packageInstalled, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dragAreaBy, hoverArea, noBalloons, noErrors, pickFromAreaContextMenu, pointerAway, wheelOverAreaTimesHolding} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridShowsColumns, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Sequences", () => {
  const session = feature(test, "features/howto/09-sequences/09-number-antibody-chains.feature", import.meta.url);
  test("Number antibody chains", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And the \"HMM\" package is installed", () => packageInstalled(page, "HMM"));
    await session.step(18, "And user opens antibodies dataset", () => openDataset(page, ds("antibodies")));
    await session.step(19, "And the Bio package is initialized", () => bioReady(page));
    await session.step(20, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(21, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(22, "And the grid shows only the columns \"Antibody_ID, AntibodyHC\"", () => gridShowsColumns(page, "Antibody_ID, AntibodyHC"));
    await session.step(23, "And the \"AntibodyHC\" column of the grid is 600 px wide", () => gridColumnWidth(page, "AntibodyHC", 600));
    await session.step(24, "When the video says \"493 antibodies, their heavy chains here\"", () => videoSays(page, "493 antibodies, their heavy chains here"));
    await session.step(26, "And user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."), undefined, "Bio > Annotate > Apply Numbering Scheme");
    await session.step(27, "Then \"Apply Antibody Numbering\" dialog should be visible", () => shouldBe(page, el("\"Apply Antibody Numbering\" dialog"), "visible"));
    await session.step(28, "And Scheme input in \"Apply Antibody Numbering\" dialog should have value \"imgt\"", () => shouldHaveValue(page, el("Scheme input in \"Apply Antibody Numbering\" dialog"), "imgt"));
    await session.step(30, "When user selects \"imgt\" in Scheme input in \"Apply Antibody Numbering\" dialog", () => selectIn(page, "imgt", el("Scheme input in \"Apply Antibody Numbering\" dialog")), undefined, "Pick a scheme: IMGT, Kabat, Chothia and more");
    await session.step(31, "And user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
    await session.step(32, "Then a new column \"AntibodyHC (aligned)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)"));
    await session.step(33, "And \"AntibodyHC (aligned)\" column should have tag \".numberingScheme\" equal to \"imgt\"", () => columnTag(page, "AntibodyHC (aligned)", ".numberingScheme", "imgt"));
    await session.step(34, "When the video says \"Frameworks and CDRs, colored on every chain\"", () => videoSays(page, "Frameworks and CDRs, colored on every chain"));
    await session.step(35, "And user hovers over monomer 30 of the \"cell 3 of AntibodyHC\" area of grid", () => hoverMonomer(page, 30, "cell 3 of AntibodyHC", el("grid")));
    await session.step(36, "Then tooltip should contain text \"CDR1\"", () => shouldContainText(page, el("tooltip"), "CDR1"));
    await session.step(37, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(38, "And the video says \"The new column: the chains aligned in IMGT positions\"", () => videoSays(page, "The new column: the chains aligned in IMGT positions"));
    await session.step(40, "And user scrolls the mouse wheel down 43 times over the \"header AntibodyHC (aligned)\" area of grid holding Shift", () => wheelOverAreaTimesHolding(page, "down", 43, "header AntibodyHC (aligned)", el("grid"), "Shift"), undefined, "Shift + wheel over its header scrolls the positions");
    await session.step(41, "And user hovers over the region track under the \"position 111A of AntibodyHC (aligned) header\" area of grid", () => hoverRegionBand(page, "position 111A of AntibodyHC (aligned) header", el("grid")));
    await session.step(42, "Then tooltip should contain text \"CDR3\"", () => shouldContainText(page, el("tooltip"), "CDR3"));
    await session.step(44, "When user picks \"Annotations > Extract CDR3 as Column\" from the context menu of the \"cell 12 of AntibodyHC (aligned)\" area of grid", () => pickFromAreaContextMenu(page, "Annotations > Extract CDR3 as Column", "cell 12 of AntibodyHC (aligned)", el("grid")), undefined, "Right-click a chain: extract CDR3 as a column");
    await session.step(45, "Then a new column \"AntibodyHC (aligned)(CDR3)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)(CDR3)"));
    await session.step(47, "When user drags the \"column resizer AntibodyHC (aligned)(CDR3)\" area of grid by 200 pixels to the right", () => dragAreaBy(page, "column resizer AntibodyHC (aligned)(CDR3)", el("grid"), 200, "right"), undefined, "Drag its edge to widen it");
    await session.step(48, "And the video says \"CDR3 alone, with its IMGT positions: 105 to 117\"", () => videoSays(page, "CDR3 alone, with its IMGT positions: 105 to 117"));
    await session.step(49, "And user hovers over the \"position 111A of AntibodyHC (aligned)(CDR3) header\" area of grid", () => hoverArea(page, "position 111A of AntibodyHC (aligned)(CDR3) header", el("grid")));
    await session.step(50, "Then tooltip should contain text \"Position: 111A\"", () => shouldContainText(page, el("tooltip"), "Position: 111A"));
    await session.step(51, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(52, "And no errors should have been logged", () => noErrors(page));
  });
});
