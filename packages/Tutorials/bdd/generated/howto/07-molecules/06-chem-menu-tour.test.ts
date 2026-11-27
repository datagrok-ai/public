/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/07-molecules/06-chem-menu-tour.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
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
import {restOnTopMenuItem} from '../../../bindings/howto-data-2.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, openDataset, packageInstalled, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Molecules", () => {
  const session = feature(test, "features/howto/07-molecules/06-chem-menu-tour.feature", import.meta.url);
  test("What the Chem menu does", {tag: ["@guide", "@tutorial", "@viewport:1280x720"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await session.step(20, "And the grid shows only the columns \"molregno, canonical_smiles, RingCount, NumHDonors, NumHAcceptors\"", () => gridShowsColumns(page, "molregno, canonical_smiles, RingCount, NumHDonors, NumHAcceptors"));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(23, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(24, "When the video says \"Chem: the tools for a table of molecules\" at the top", () => videoSaysAt(page, "Chem: the tools for a table of molecules", "top"));
    await session.step(26, "And user rests the pointer on \"Chem > Search\" in the top menu", () => restOnTopMenuItem(page, "Chem > Search"), undefined, "Search: by a fragment, by similarity, or for diversity");
    await session.step(27, "And user rests the pointer on \"Chem > Search > Similarity Search...\" in the top menu", () => restOnTopMenuItem(page, "Chem > Search > Similarity Search..."));
    await session.step(28, "And user rests the pointer on \"Chem > Search > Substructure Search...\" in the top menu", () => restOnTopMenuItem(page, "Chem > Search > Substructure Search..."));
    await session.step(29, "Then \"Substructure Search...\" menu item should be enabled", () => shouldBe(page, el("\"Substructure Search...\" menu item"), "enabled"));
    await session.step(31, "When user rests the pointer on \"Chem > Calculate\" in the top menu", () => restOnTopMenuItem(page, "Chem > Calculate"), undefined, "Calculate: properties, scores and alerts as new columns");
    await session.step(32, "And user rests the pointer on \"Chem > Calculate > Toxicity Risks...\" in the top menu", () => restOnTopMenuItem(page, "Chem > Calculate > Toxicity Risks..."));
    await session.step(33, "And user rests the pointer on \"Chem > Calculate > MPO Score...\" in the top menu", () => restOnTopMenuItem(page, "Chem > Calculate > MPO Score..."));
    await session.step(34, "And user rests the pointer on \"Chem > Calculate > Chemical Properties...\" in the top menu", () => restOnTopMenuItem(page, "Chem > Calculate > Chemical Properties..."));
    await session.step(35, "Then \"Chemical Properties...\" menu item should be enabled", () => shouldBe(page, el("\"Chemical Properties...\" menu item"), "enabled"));
    await session.step(37, "When user rests the pointer on \"Chem > Analyze\" in the top menu", () => restOnTopMenuItem(page, "Chem > Analyze"), undefined, "Analyze: how structure relates to activity");
    await session.step(38, "And user rests the pointer on \"Chem > Analyze > Chemical Space...\" in the top menu", () => restOnTopMenuItem(page, "Chem > Analyze > Chemical Space..."));
    await session.step(39, "And user rests the pointer on \"Chem > Analyze > R-Groups Analysis...\" in the top menu", () => restOnTopMenuItem(page, "Chem > Analyze > R-Groups Analysis..."));
    await session.step(40, "And user rests the pointer on \"Chem > Analyze > Activity Cliffs...\" in the top menu", () => restOnTopMenuItem(page, "Chem > Analyze > Activity Cliffs..."));
    await session.step(41, "And user rests the pointer on \"Chem > Analyze > Scaffold Tree\" in the top menu", () => restOnTopMenuItem(page, "Chem > Analyze > Scaffold Tree"));
    await session.step(42, "Then \"Activity Cliffs...\" menu item should be enabled", () => shouldBe(page, el("\"Activity Cliffs...\" menu item"), "enabled"));
    await session.step(44, "When user rests the pointer on \"Chem > Transform\" in the top menu", () => restOnTopMenuItem(page, "Chem > Transform"), undefined, "Transform: change the molecules: notation, cleanup, reactions");
    await session.step(45, "And user rests the pointer on \"Chem > Transform > Convert Notation...\" in the top menu", () => restOnTopMenuItem(page, "Chem > Transform > Convert Notation..."));
    await session.step(46, "And user rests the pointer on \"Chem > Transform > Reactions\" in the top menu", () => restOnTopMenuItem(page, "Chem > Transform > Reactions"));
    await session.step(47, "Then \"Convert Notation...\" menu item should be enabled", () => shouldBe(page, el("\"Convert Notation...\" menu item"), "enabled"));
    await session.step(49, "When user rests the pointer on \"Chem > Admetica\" in the top menu", () => restOnTopMenuItem(page, "Chem > Admetica"), undefined, "Admetica: ADMET predictions, computed on the server");
    await session.step(50, "Then \"Admetica\" menu item should be enabled", () => shouldBe(page, el("\"Admetica\" menu item"), "enabled"));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
});
