/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/chem/activity-cliffs.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.activity-cliffs]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {startTutorial, stepDone, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, toggle} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, noHintShown, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickFirstArea, clickTakenArea, dragZoomOverArea, hoverFirstArea, hoverFirstFreeArea, noErrors, readingIs, readingReads, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragTopBorder} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Activity Cliffs tutorial", () => {
  const session = feature(test, "features/chem/activity-cliffs.feature", import.meta.url);
  test("A learner completes the Activity Cliffs tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.activity-cliffs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(23, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(24, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(25, "And the \"Activity Cliffs\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Activity Cliffs"));
    await session.step(26, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(29, "When user starts the \"Activity Cliffs\" tutorial", () => startTutorial(page, "Activity Cliffs"));
    await session.step(30, "Then the tutorial progress should be 1 of 12", () => tutorialProgress(page, 1, 12));
    await session.step(31, "When user picks \"Chem > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Activity Cliffs..."));
    await session.step(32, "Then the tutorial step \"On the Top Menu, click Chem > Analyze > Activity Cliffs...\" should be done", () => stepDone(page, "On the Top Menu, click Chem > Analyze > Activity Cliffs..."));
    await session.step(33, "And \"Activity Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Activity Cliffs\" dialog"), "visible"));
    await session.step(34, "When user clicks on OK button in \"Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Activity Cliffs\" dialog")));
    await session.step(35, "Then the tutorial step \"Click OK\" should be done", () => stepDone(page, "Click OK"));
    await session.step(36, "And the tutorial step \"Wait for analysis to complete\" should be done", () => stepDone(page, "Wait for analysis to complete"));
    await session.step(37, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(38, "And the \"cliffs\" reading of scatter plot viewer should be 15", () => readingIs(page, "cliffs", el("scatter plot viewer"), 15));
    await session.step(40, "When user hovers over the first \"marker\" area of scatter plot viewer", () => hoverFirstArea(page, "marker", el("scatter plot viewer")));
    await session.step(41, "Then the tutorial step \"Hover over data points for molecule information\" should be done", () => stepDone(page, "Hover over data points for molecule information"));
    await session.step(42, "And tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
    await session.step(43, "When user toggles \"Show only cliffs\" input in scatter plot viewer", () => toggle(page, el("\"Show only cliffs\" input in scatter plot viewer")));
    await session.step(44, "Then the tutorial step \"To view only the cliffs, toggle Show only cliffs.\" should be done", () => stepDone(page, "To view only the cliffs, toggle Show only cliffs."));
    await session.step(45, "And the \"only cliffs\" reading of scatter plot viewer should be \"true\"", () => readingReads(page, "only cliffs", el("scatter plot viewer"), "true"));
    await session.step(47, "When user drags a zoom box over the \"view\" area of scatter plot viewer", () => dragZoomOverArea(page, "view", el("scatter plot viewer")));
    await session.step(48, "Then the tutorial step \"Press Use Alt + Mouse Drag to zoom in\" should be done", () => stepDone(page, "Press Use Alt + Mouse Drag to zoom in"));
    await session.step(50, "When user hovers over the first free \"line\" area of scatter plot viewer", () => hoverFirstFreeArea(page, "line", el("scatter plot viewer")));
    await session.step(51, "Then the tutorial step \"Hover over the green line to see the pair of molecules\" should be done", () => stepDone(page, "Hover over the green line to see the pair of molecules"));
    await session.step(52, "And tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
    await session.step(53, "When user clicks on that area of scatter plot viewer", () => clickTakenArea(page, el("scatter plot viewer")));
    await session.step(54, "Then the tutorial step \"Click on the green line connecting that molecule pair\" should be done", () => stepDone(page, "Click on the green line connecting that molecule pair"));
    await session.step(55, "And \"Cliff Details\" pane in context panel should be visible", () => shouldBe(page, el("\"Cliff Details\" pane in context panel"), "visible"));
    await session.step(58, "When user clicks on second cliff molecule in context panel", () => clickOn(page, el("second cliff molecule in context panel")));
    await session.step(59, "Then the tutorial step \"On the Context Panel, click any molecule\" should be done", () => stepDone(page, "On the Context Panel, click any molecule"));
    await session.step(61, "When user clicks on \"15 cliffs\" button in scatter plot viewer", () => clickOn(page, el("\"15 cliffs\" button in scatter plot viewer")));
    await session.step(62, "Then the tutorial step \"At the top right corner of the scatterplot, click 15 CLIFFS\" should be done", () => stepDone(page, "At the top right corner of the scatterplot, click 15 CLIFFS"));
    await session.step(63, "And the open tableview should have 2 grid viewers", () => viewerCount(page, 2, "grid"));
    await session.step(64, "When user drags the top border of second grid viewer by 150 pixels up", () => dragTopBorder(page, el("second grid viewer"), 150));
    await session.step(65, "Then the tutorial step \"Drag the top border of that table upwards to create more space for it\" should be done", () => stepDone(page, "Drag the top border of that table upwards to create more space for it"));
    await session.step(66, "When user clicks on the first \"cell\" area of second grid viewer", () => clickFirstArea(page, "cell", el("second grid viewer")));
    await session.step(67, "Then the tutorial step \"In the cliffs table, click any cell in the first row\" should be done", () => stepDone(page, "In the cliffs table, click any cell in the first row"));
    await session.step(69, "And the \"Activity Cliffs\" tutorial should be completed", () => tutorialCompleted(page, "Activity Cliffs"));
    await session.step(70, "And the tutorial should have listed 12 steps", () => tutorialStepsListed(page, 12));
    await session.step(71, "And the tutorial progress should be 12 of 12", () => tutorialProgress(page, 12, 12));
    await session.step(72, "And no hint should be shown", () => noHintShown(page));
    await session.step(73, "And no errors should have been logged", () => noErrors(page));
  });
});
