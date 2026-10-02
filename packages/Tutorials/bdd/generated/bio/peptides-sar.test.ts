/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/bio/peptides-sar.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.peptides-sar]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {startTutorial, stepDone, stepDoneTimes, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, enterInto, hoverOver, pressKey, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {noneSelected, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {customFired, listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {autostartsCompleted, contextPanelShows, noHintShown, packageInstalled, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverArea, noErrors, propertyShouldBe, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Peptides SAR tutorial", () => {
  const session = feature(test, "features/bio/peptides-sar.feature", import.meta.url);
  test("A learner completes the Peptides SAR tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.peptides-sar"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the \"Bio\" package is installed", () => packageInstalled(page, "Bio"));
    await session.step(23, "And the \"Peptides\" package is installed", () => packageInstalled(page, "Peptides"));
    await session.step(24, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(25, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(26, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(27, "And the \"Peptides SAR\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Peptides SAR"));
    await session.step(28, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(31, "When user starts the \"Peptides SAR\" tutorial", () => startTutorial(page, "Peptides SAR"));
    await session.step(32, "Then the tutorial progress should be 1 of 24", () => tutorialProgress(page, 1, 24));
    await session.step(33, "When user picks \"Bio > Analyze > SAR...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > SAR..."));
    await session.step(34, "Then the tutorial step \"On the Top Menu, click Bio > Analyze > SAR...\" should be done", () => stepDone(page, "On the Top Menu, click Bio > Analyze > SAR..."));
    await session.step(35, "And \"Analyze Peptides\" dialog should be visible", () => shouldBe(page, el("\"Analyze Peptides\" dialog"), "visible"));
    await session.step(36, "When user clicks on \"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog", () => clickOn(page, el("\"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog")));
    await session.step(37, "Then the tutorial step \"Click the Gear Icon (⚙)\" should be done", () => stepDone(page, "Click the Gear Icon (⚙)"));
    await session.step(38, "When user enters \"90\" into \"Similarity Threshold\" input in \"Analyze Peptides\" dialog", () => enterInto(page, "90", el("\"Similarity Threshold\" input in \"Analyze Peptides\" dialog")));
    await session.step(39, "Then the tutorial step \"Set Similarity Threshold to 90\" should be done", () => stepDone(page, "Set Similarity Threshold to 90"));
    await session.step(40, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(41, "When user clicks on OK button in \"Analyze Peptides\" dialog", () => clickOn(page, el("OK button in \"Analyze Peptides\" dialog")));
    await session.step(42, "Then the tutorial step \"Click OK to start analysis\" should be done", () => stepDone(page, "Click OK to start analysis"));
    await session.step(43, "And the \"peptides-sar-ready\" custom event should have fired", () => customFired(page, "peptides-sar-ready"));
    await session.step(44, "And the tutorial step \"Wait for analysis to complete\" should be done", () => stepDone(page, "Wait for analysis to complete"));
    await session.step(45, "When user clicks on NEXT button in hint popup", () => clickOn(page, el("NEXT button in hint popup")));
    await session.step(46, "Then the tutorial step \"Click NEXT to proceed\" should be done", () => stepDone(page, "Click NEXT to proceed"));
    await session.step(48, "When user hovers over the \"cell 1 of 2\" area of grid", () => hoverArea(page, "cell 1 of 2", el("grid")));
    await session.step(49, "Then the tutorial step \"Hover over a monomer cell in the main table grid\" should be done", () => stepDone(page, "Hover over a monomer cell in the main table grid"));
    await session.step(50, "And tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
    await session.step(51, "When user clicks on the \"Aca at 3\" area of grid", () => clickArea(page, "Aca at 3", el("grid")));
    await session.step(52, "Then the tutorial step \"Click any WebLogo letter\" should be done", () => stepDone(page, "Click any WebLogo letter"));
    await session.step(53, "And some rows should be selected", () => someSelected(page));
    await session.step(54, "When user clicks on NEXT button in hint popup", () => clickOn(page, el("NEXT button in hint popup")));
    await session.step(55, "Then the tutorial step \"Click NEXT to proceed\" should be done 2 times", () => stepDoneTimes(page, "Click NEXT to proceed", 2));
    await session.step(56, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(57, "Then the tutorial step \"Press Esc to clear selection\" should be done", () => stepDone(page, "Press Esc to clear selection"));
    await session.step(58, "And no rows should be selected", () => noneSelected(page));
    await session.step(60, "When user hovers over Sequence Variability Map viewer", () => hoverOver(page, el("Sequence Variability Map viewer")));
    await session.step(61, "And user clicks on settings icon of Sequence Variability Map viewer", () => clickOn(page, el("settings icon of Sequence Variability Map viewer")));
    await session.step(62, "Then the tutorial step \"Open SVM settings (gear)\" should be done", () => stepDone(page, "Open SVM settings (gear)"));
    await session.step(63, "And the context panel should show \"Sequence Variability Map\"", () => contextPanelShows(page, "Sequence Variability Map"));
    await session.step(64, "When user clicks on NEXT button in hint popup", () => clickOn(page, el("NEXT button in hint popup")));
    await session.step(65, "Then the tutorial step \"Click NEXT to proceed\" should be done 3 times", () => stepDoneTimes(page, "Click NEXT to proceed", 3));
    await session.step(67, "When user clicks on the \"cell 1Nal at 1\" area of Sequence Variability Map viewer", () => clickArea(page, "cell 1Nal at 1", el("Sequence Variability Map viewer")));
    await session.step(68, "Then the tutorial step \"Click a Mutation Cliffs cell\" should be done", () => stepDone(page, "Click a Mutation Cliffs cell"));
    await session.step(69, "And \"Mutation Cliffs pairs\" pane in context panel should be visible", () => shouldBe(page, el("\"Mutation Cliffs pairs\" pane in context panel"), "visible"));
    await session.step(70, "When user clicks on NEXT button in hint popup", () => clickOn(page, el("NEXT button in hint popup")));
    await session.step(71, "Then the tutorial step \"Click NEXT to proceed\" should be done 4 times", () => stepDoneTimes(page, "Click NEXT to proceed", 4));
    await session.step(72, "When user clicks on \"Invariant Map\" checkbox in Sequence Variability Map viewer", () => clickOn(page, el("\"Invariant Map\" checkbox in Sequence Variability Map viewer")));
    await session.step(73, "Then the tutorial step \"Switch SVM mode to Invariant Map\" should be done", () => stepDone(page, "Switch SVM mode to Invariant Map"));
    await session.step(74, "And the \"mode\" reading of Sequence Variability Map viewer should be \"Invariant Map\"", () => readingReads(page, "mode", el("Sequence Variability Map viewer"), "Invariant Map"));
    await session.step(75, "When user clicks on NEXT button in hint popup", () => clickOn(page, el("NEXT button in hint popup")));
    await session.step(76, "Then the tutorial step \"Click NEXT to proceed\" should be done 5 times", () => stepDoneTimes(page, "Click NEXT to proceed", 5));
    await session.step(78, "When user hovers over Most Potent Residues viewer", () => hoverOver(page, el("Most Potent Residues viewer")));
    await session.step(79, "And user clicks on settings icon of Most Potent Residues viewer", () => clickOn(page, el("settings icon of Most Potent Residues viewer")));
    await session.step(80, "Then the tutorial step \"Open Most Potent Residues settings (gear)\" should be done", () => stepDone(page, "Open Most Potent Residues settings (gear)"));
    await session.step(81, "And the context panel should show \"Most Potent Residues\"", () => contextPanelShows(page, "Most Potent Residues"));
    await session.step(82, "When user clicks on NEXT button in hint popup", () => clickOn(page, el("NEXT button in hint popup")));
    await session.step(83, "Then the tutorial step \"Click NEXT to proceed\" should be done 6 times", () => stepDoneTimes(page, "Click NEXT to proceed", 6));
    await session.step(84, "When user hovers over Logo Summary Table viewer", () => hoverOver(page, el("Logo Summary Table viewer")));
    await session.step(85, "And user clicks on settings icon of Logo Summary Table viewer", () => clickOn(page, el("settings icon of Logo Summary Table viewer")));
    await session.step(86, "Then the tutorial step \"Open Logo Summary Table settings (gear)\" should be done", () => stepDone(page, "Open Logo Summary Table settings (gear)"));
    await session.step(87, "And the context panel should show \"Logo Summary Table\"", () => contextPanelShows(page, "Logo Summary Table"));
    await session.step(88, "When user clicks on \"...\" button in \"Columns\" property", () => clickOn(page, el("\"...\" button in \"Columns\" property")));
    await session.step(89, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(90, "When user types \"14\" into \"Search\" input in \"Select columns...\" dialog", () => typeInto(page, "14", el("\"Search\" input in \"Select columns...\" dialog")));
    await session.step(91, "And user toggles the \"14\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "14", el("\"Select columns...\" dialog")));
    await session.step(92, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(93, "Then the tutorial step \"Add pie chart aggregation for position 14\" should be done", () => stepDone(page, "Add pie chart aggregation for position 14"));
    await session.step(94, "And \"columns\" property of Logo Summary Table viewer should be \"14\"", () => propertyShouldBe(page, "columns", el("Logo Summary Table viewer"), "14"));
    await session.step(95, "When user clicks on NEXT button in hint popup", () => clickOn(page, el("NEXT button in hint popup")));
    await session.step(96, "Then the tutorial step \"Scroll horizontally in Logo Summary Table. Click NEXT to proceed to next step\" should be done", () => stepDone(page, "Scroll horizontally in Logo Summary Table. Click NEXT to proceed to next step"));
    await session.step(98, "When user clicks on \"Peptides analysis settings\" icon", () => clickOn(page, el("\"Peptides analysis settings\" icon")));
    await session.step(99, "Then the tutorial step \"Click the Wrench to update analysis configuration\" should be done", () => stepDone(page, "Click the Wrench to update analysis configuration"));
    await session.step(100, "And \"Peptides settings\" dialog should be visible", () => shouldBe(page, el("\"Peptides settings\" dialog"), "visible"));
    await session.step(101, "When user checks \"Dendrogram\" input in \"Peptides settings\" dialog", () => check(page, el("\"Dendrogram\" input in \"Peptides settings\" dialog")));
    await session.step(102, "Then the tutorial step \"Check Dendrogram\" should be done", () => stepDone(page, "Check Dendrogram"));
    await session.step(103, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(104, "When user clicks on OK button in \"Peptides settings\" dialog", () => clickOn(page, el("OK button in \"Peptides settings\" dialog")));
    await session.step(105, "Then the tutorial step \"Click OK to re-run analysis\" should be done", () => stepDone(page, "Click OK to re-run analysis"));
    await session.step(106, "And the \"peptides-sar-ready\" custom event should have fired", () => customFired(page, "peptides-sar-ready"));
    await session.step(107, "When user clicks on OK button in hint popup", () => clickOn(page, el("OK button in hint popup")));
    await session.step(109, "And the \"Peptides SAR\" tutorial should be completed", () => tutorialCompleted(page, "Peptides SAR"));
    await session.step(110, "And the tutorial should have listed 24 steps", () => tutorialStepsListed(page, 24));
    await session.step(111, "And the tutorial progress should be 24 of 24", () => tutorialProgress(page, 24, 24));
    await session.step(112, "And no hint should be shown", () => noHintShown(page));
    await session.step(113, "And no errors should have been logged", () => noErrors(page));
  });
});
