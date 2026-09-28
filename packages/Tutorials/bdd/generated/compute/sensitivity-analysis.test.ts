/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/compute/sensitivity-analysis.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.sensitivity-analysis]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {startTutorial, stepDone, stepDoneTimes, stepListedTimes, stepNotDone, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen, walkTour} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, selectIn, switchOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, noHintShown, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dragAreaToArea, noErrors, readingIs, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Sensitivity Analysis tutorial", () => {
  const session = feature(test, "features/compute/sensitivity-analysis.feature", import.meta.url);
  test("A learner completes the Sensitivity analysis tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.sensitivity-analysis"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(21, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(22, "And the \"Sensitivity analysis\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Sensitivity analysis"));
    await session.step(23, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(26, "When user starts the \"Sensitivity analysis\" tutorial", () => startTutorial(page, "Sensitivity analysis"));
    await session.step(27, "Then the tutorial progress should be 1 of 16", () => tutorialProgress(page, 1, 16));
    await session.step(28, "When user clicks on Apps tree node inside browse tree", () => clickOn(page, el("Apps tree node inside browse tree")));
    await session.step(29, "Then the tutorial step \"Open Apps\" should be done", () => stepDone(page, "Open Apps"));
    await session.step(30, "Given the tutorial step \"Run Model Hub\" should not be done yet", () => stepNotDone(page, "Run Model Hub"));
    await session.step(31, "When user double-clicks on Model-Hub gallery card", () => doubleClickOn(page, el("Model-Hub gallery card")));
    await session.step(32, "Then the tutorial step \"Run Model Hub\" should be done", () => stepDone(page, "Run Model Hub"));
    await session.step(33, "Given the tutorial step \"Run the \\\"Ball flight\\\" model\" should not be done yet", () => stepNotDone(page, "Run the \"Ball flight\" model"));
    await session.step(34, "When user double-clicks on \"Ball Flight Simulation\" link", () => doubleClickOn(page, el("\"Ball Flight Simulation\" link")));
    await session.step(35, "Then the tutorial step \"Run the \\\"Ball flight\\\" model\" should be done", () => stepDone(page, "Run the \"Ball flight\" model"));
    await session.step(36, "Given the tutorial step \"Click \\\"OK\\\"\" should not be done yet", () => stepNotDone(page, "Click \"OK\""));
    await session.step(37, "When user goes through the tour to its end", () => walkTour(page));
    await session.step(38, "Then the tutorial step \"Click \\\"OK\\\"\" should be done", () => stepDone(page, "Click \"OK\""));
    await session.step(40, "When user clicks on \"Run sensitivity analysis\" icon", () => clickOn(page, el("\"Run sensitivity analysis\" icon")));
    await session.step(41, "Then the tutorial step \"Run sensitivity analysis\" should be done", () => stepDone(page, "Run sensitivity analysis"));
    await session.step(42, "Given the tutorial step \"Set \\\"Samples\\\" to 100\" should not be done yet", () => stepNotDone(page, "Set \"Samples\" to 100"));
    await session.step(43, "When user enters \"100\" into \"Samples\" input", () => enterInto(page, "100", el("\"Samples\" input")));
    await session.step(44, "Then the tutorial step \"Set \\\"Samples\\\" to 100\" should be done", () => stepDone(page, "Set \"Samples\" to 100"));
    await session.step(45, "When user switches on \"Angle\" input", () => switchOn(page, el("\"Angle\" input")));
    await session.step(46, "Then the tutorial step \"Toggle the \\\"Angle\\\" parameter\" should be done", () => stepDone(page, "Toggle the \"Angle\" parameter"));
    await session.step(47, "When user clicks on \"Run\" icon", () => clickOn(page, el("\"Run\" icon")));
    await session.step(48, "Then the tutorial step \"Run sensitivity analysis\" should be done 2 times", () => stepDoneTimes(page, "Run sensitivity analysis", 2));
    await session.step(49, "Given the tutorial step \"Explore each viewer\" should not be done yet", () => stepNotDone(page, "Explore each viewer"));
    await session.step(50, "When user goes through the tour to its end", () => walkTour(page));
    await session.step(51, "Then the tutorial step \"Explore each viewer\" should be done", () => stepDone(page, "Explore each viewer"));
    await session.step(53, "When user drags the \"range min handle \\\"maxDist\\\"\" area of pc plot viewer to the \"range max handle \\\"maxDist\\\"\" area", () => dragAreaToArea(page, "range min handle \"maxDist\"", el("pc plot viewer"), "range max handle \"maxDist\""));
    await session.step(54, "Then the tutorial step \"Move slider\" should be done", () => stepDone(page, "Move slider"));
    await session.step(56, "And the \"rows shown\" reading of pc plot viewer should be 1", () => readingIs(page, "rows shown", el("pc plot viewer"), 1));
    await session.step(57, "Given the tutorial step \"Explore the solution\" should not be done yet", () => stepNotDone(page, "Explore the solution"));
    await session.step(58, "When user goes through the tour to its end", () => walkTour(page));
    await session.step(59, "Then the tutorial step \"Explore the solution\" should be done", () => stepDone(page, "Explore the solution"));
    await session.step(61, "When user selects \"Sobol\" in \"Method\" input", () => selectIn(page, "Sobol", el("\"Method\" input")));
    await session.step(62, "Then the tutorial step \"Set \\\"Method\\\" to \\\"Sobol\\\"\" should be done", () => stepDone(page, "Set \"Method\" to \"Sobol\""));
    await session.step(63, "When user switches on \"Velocity\" input", () => switchOn(page, el("\"Velocity\" input")));
    await session.step(64, "Then the tutorial step \"Toggle the \\\"Velocity\\\" parameter\" should be done", () => stepDone(page, "Toggle the \"Velocity\" parameter"));
    await session.step(65, "When user clicks on \"Run\" icon", () => clickOn(page, el("\"Run\" icon")));
    await session.step(66, "Then the tutorial step \"Run sensitivity analysis\" should be done 3 times", () => stepDoneTimes(page, "Run sensitivity analysis", 3));
    await session.step(67, "And the open tableview should have 2 bar chart viewers", () => viewerCount(page, 2, "bar chart"));
    await session.step(68, "Given the tutorial step \"Explore each viewer\" should be listed 2 times", () => stepListedTimes(page, "Explore each viewer", 2));
    await session.step(69, "When user goes through the tour to its end", () => walkTour(page));
    await session.step(70, "Then the tutorial step \"Explore each viewer\" should be done 2 times", () => stepDoneTimes(page, "Explore each viewer", 2));
    await session.step(71, "Given the tutorial step \"Click \\\"Clear\\\"\" should not be done yet", () => stepNotDone(page, "Click \"Clear\""));
    await session.step(72, "When user goes through the tour to its end", () => walkTour(page));
    await session.step(73, "Then the tutorial step \"Click \\\"Clear\\\"\" should be done", () => stepDone(page, "Click \"Clear\""));
    await session.step(75, "And the \"Sensitivity analysis\" tutorial should be completed", () => tutorialCompleted(page, "Sensitivity analysis"));
    await session.step(76, "And the tutorial should have listed 16 steps", () => tutorialStepsListed(page, 16));
    await session.step(77, "And the tutorial progress should be 16 of 16", () => tutorialProgress(page, 16, 16));
    await session.step(78, "And no hint should be shown", () => noHintShown(page));
    await session.step(79, "And no errors should have been logged", () => noErrors(page));
  });
});
