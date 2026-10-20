/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/chem/similarity-diversity.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.similarity-diversity-search]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {startTutorial, stepDone, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, hoverOver, pressKeyIn, shouldBe, typeInto, uncheck, waitSeconds} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {colorCodedAs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, contextPanelShows, contextPanelShowsCurrentCell, noHintShown, packageInstalled, sketcherIs, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, pickFromAreaContextMenu, pickFromOpenMenu, propertyShouldBe, readingNotAsRemembered, readingReads, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {clickNthCardBesideDrawing, clickOtherCell, currentRowIsClickedCard, hoverNthCard, toggleInColumnList, walkToColumn} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Similarity and Diversity Search tutorial", () => {
  const session = feature(test, "features/chem/similarity-diversity.feature", import.meta.url);
  test("A learner completes the Similarity and Diversity Search tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.similarity-diversity-search"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(31, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(32, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(33, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(34, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(35, "And the \"Similarity and Diversity Search\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Similarity and Diversity Search"));
    await session.step(36, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(39, "When user starts the \"Similarity and Diversity Search\" tutorial", () => startTutorial(page, "Similarity and Diversity Search"));
    await session.step(40, "Then the tutorial progress should be 1 of 12", () => tutorialProgress(page, 1, 12));
    await session.step(41, "When user picks \"Chem > Search > Similarity Search...\" from the top menu", () => pickFromTopMenu(page, "Chem > Search > Similarity Search..."));
    await session.step(42, "Then the tutorial step \"On the Top Menu, click Chem > Search > Similarity Search...\" should be done", () => stepDone(page, "On the Top Menu, click Chem > Search > Similarity Search..."));
    await session.step(43, "And Chem Similarity Search viewer should be visible", () => shouldBe(page, el("Chem Similarity Search viewer"), "visible"));
    await session.step(44, "When user picks \"Chem > Search > Diversity Search...\" from the top menu", () => pickFromTopMenu(page, "Chem > Search > Diversity Search..."));
    await session.step(45, "Then the tutorial step \"Next, click Chem > Search > Diversity Search...\" should be done", () => stepDone(page, "Next, click Chem > Search > Diversity Search..."));
    await session.step(46, "And Chem Diversity Search viewer should be visible", () => shouldBe(page, el("Chem Diversity Search viewer"), "visible"));
    await session.step(48, "When user clicks on card 2 of Chem Similarity Search viewer beside its drawing", () => clickNthCardBesideDrawing(page, 2, el("Chem Similarity Search viewer")));
    await session.step(49, "Then the tutorial step \"On the Most similar structures viewer, click the molecule next to the reference molecule\" should be done", () => stepDone(page, "On the Most similar structures viewer, click the molecule next to the reference molecule"));
    await session.step(50, "And the current row should be the row of the clicked card", () => currentRowIsClickedCard(page));
    await session.step(51, "When user clicks on card 3 of Chem Diversity Search viewer beside its drawing", () => clickNthCardBesideDrawing(page, 3, el("Chem Diversity Search viewer")));
    await session.step(52, "Then the tutorial step \"Now, click any molecule in the diversity viewer\" should be done", () => stepDone(page, "Now, click any molecule in the diversity viewer"));
    await session.step(53, "And the current row should be the row of the clicked card", () => currentRowIsClickedCard(page));
    await session.step(56, "When user waits 1 second", () => waitSeconds(page, 1));
    await session.step(58, "When user hovers over Chem Similarity Search viewer", () => hoverOver(page, el("Chem Similarity Search viewer")));
    await session.step(59, "And user clicks on settings icon of Chem Similarity Search viewer", () => clickOn(page, el("settings icon of Chem Similarity Search viewer")));
    await session.step(60, "Then the tutorial step \"Hover over similarity viewer and click gear icon in the right top corner of the viewer to open settings\" should be done", () => stepDone(page, "Hover over similarity viewer and click gear icon in the right top corner of the viewer to open settings"));
    await session.step(61, "And the context panel should show \"Chem Similarity Search\"", () => contextPanelShows(page, "Chem Similarity Search"));
    await session.step(62, "When user expands Misc category", () => expand(page, el("Misc category")));
    await session.step(63, "And user unchecks \"Follow Current Row\" property", () => uncheck(page, el("\"Follow Current Row\" property")));
    await session.step(64, "Then the tutorial step \"Under Misc, clear the Follow Current Row checkbox\" should be done", () => stepDone(page, "Under Misc, clear the Follow Current Row checkbox"));
    await session.step(65, "And \"followCurrentRow\" property of Chem Similarity Search viewer should be \"false\"", () => propertyShouldBe(page, "followCurrentRow", el("Chem Similarity Search viewer"), "false"));
    await session.step(67, "When user remembers the \"scores\" reading of Chem Similarity Search viewer", () => rememberReading(page, "scores", el("Chem Similarity Search viewer")));
    await session.step(68, "And user clicks on \"Edit\" icon in Chem Similarity Search viewer", () => clickOn(page, el("\"Edit\" icon in Chem Similarity Search viewer")));
    await session.step(69, "Then the tutorial step \"On the reference molecule, click the Edit icon\" should be done", () => stepDone(page, "On the reference molecule, click the Edit icon"));
    await session.step(70, "And sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(71, "When user types \"CNc1nc(Nc2ccc(Br)cc2)nc(N)c1[N+](=O)[O-]\" into molecule input of sketcher dialog", () => typeInto(page, "CNc1nc(Nc2ccc(Br)cc2)nc(N)c1[N+](=O)[O-]", el("molecule input of sketcher dialog")));
    await session.step(72, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(73, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(74, "Then the tutorial step \"Set new reference molecule\" should be done", () => stepDone(page, "Set new reference molecule"));
    await session.step(75, "And the \"scores\" reading of Chem Similarity Search viewer should not be as remembered", () => readingNotAsRemembered(page, "scores", el("Chem Similarity Search viewer")));
    await session.step(77, "When user clicks on a \"smiles\" cell of grid other than the current one", () => clickOtherCell(page, "smiles", el("grid")));
    await session.step(78, "Then the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
    await session.step(79, "When user hovers over card 1 of Chem Similarity Search viewer", () => hoverNthCard(page, 1, el("Chem Similarity Search viewer")));
    await session.step(80, "And user clicks on \"More\" icon in Chem Similarity Search viewer", () => clickOn(page, el("\"More\" icon in Chem Similarity Search viewer")));
    await session.step(81, "And user picks \"Explore\" from the open menu", () => pickFromOpenMenu(page, "Explore"));
    await session.step(82, "Then the tutorial step \"Hover over the reference molecule, click the More icon, and then Explore\" should be done", () => stepDone(page, "Hover over the reference molecule, click the More icon, and then Explore"));
    await session.step(86, "When user waits 1 second", () => waitSeconds(page, 1));
    await session.step(87, "And user clicks on \"Tanimoto, Morgan\" link in Chem Diversity Search viewer", () => clickOn(page, el("\"Tanimoto, Morgan\" link in Chem Diversity Search viewer")));
    await session.step(88, "Then the tutorial step \"In the top right corner of the diversity viewer, click Tanimoto, Morgan\" should be done", () => stepDone(page, "In the top right corner of the diversity viewer, click Tanimoto, Morgan"));
    await session.step(89, "And the context panel should show \"Chem Diversity Search\"", () => contextPanelShows(page, "Chem Diversity Search"));
    await session.step(90, "When user clicks on \"...\" button in \"Molecule Properties\" property", () => clickOn(page, el("\"...\" button in \"Molecule Properties\" property")));
    await session.step(91, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(92, "When user types \"NumValenceElectrons\" into \"Search\" input in \"Select columns...\" dialog", () => typeInto(page, "NumValenceElectrons", el("\"Search\" input in \"Select columns...\" dialog")));
    await session.step(93, "And user toggles the \"NumValenceElectrons\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "NumValenceElectrons", el("\"Select columns...\" dialog")));
    await session.step(94, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(95, "Then the tutorial step \"Select NumValenceElectrons column\" should be done", () => stepDone(page, "Select NumValenceElectrons column"));
    await session.step(96, "And the \"card properties\" reading of Chem Diversity Search viewer should be \"NumValenceElectrons\"", () => readingReads(page, "card properties", el("Chem Diversity Search viewer"), "NumValenceElectrons"));
    await session.step(98, "When user clicks on the \"header smiles\" area of grid", () => clickArea(page, "header smiles", el("grid")));
    await session.step(99, "And user moves the current cell of grid to the \"NumValenceElectrons\" column", () => walkToColumn(page, el("grid"), "NumValenceElectrons"));
    await session.step(100, "And user picks \"Color Coding > Linear\" from the context menu of the \"header NumValenceElectrons\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header NumValenceElectrons", el("grid")));
    await session.step(101, "Then the tutorial step \"Add Color Coding for NumValenceElectrons column\" should be done", () => stepDone(page, "Add Color Coding for NumValenceElectrons column"));
    await session.step(102, "And \"NumValenceElectrons\" column should be color-coded linearly", () => colorCodedAs(page, "NumValenceElectrons", "linearly"));
    await session.step(104, "And the \"Similarity and Diversity Search\" tutorial should be completed", () => tutorialCompleted(page, "Similarity and Diversity Search"));
    await session.step(105, "And the tutorial should have listed 12 steps", () => tutorialStepsListed(page, 12));
    await session.step(106, "And the tutorial progress should be 12 of 12", () => tutorialProgress(page, 12, 12));
    await session.step(107, "And no hint should be shown", () => noHintShown(page));
    await session.step(108, "And no errors should have been logged", () => noErrors(page));
  });
});
