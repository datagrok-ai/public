/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/compute/parameter-optimization.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.parameter-optimization]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {startTutorial, stepDone, stepDoneTimes, stepNotDone, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen, walkTour} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, selectIn, shouldBe, switchOff, switchOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, noHintShown, packageInstalled, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Parameter Optimization tutorial", () => {
  const session = feature(test, "features/compute/parameter-optimization.feature", import.meta.url);
  test("A learner completes the Parameter optimization tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.parameter-optimization"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the \"Compute2\" package is installed", () => packageInstalled(page, "Compute2"));
    await session.step(16, "And the \"DiffStudio\" package is installed", () => packageInstalled(page, "DiffStudio"));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(19, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(20, "And the \"Parameter optimization\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Parameter optimization"));
    await session.step(21, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(24, "When user starts the \"Parameter optimization\" tutorial", () => startTutorial(page, "Parameter optimization"));
    await session.step(25, "Then the tutorial progress should be 1 of 15", () => tutorialProgress(page, 1, 15));
    await session.step(26, "When user clicks on Apps tree node inside browse tree", () => clickOn(page, el("Apps tree node inside browse tree")));
    await session.step(27, "Then the tutorial step \"Open Apps\" should be done", () => stepDone(page, "Open Apps"));
    await session.step(28, "Given the tutorial step \"Run Model Hub\" should not be done yet", () => stepNotDone(page, "Run Model Hub"));
    await session.step(29, "When user double-clicks on Model-Hub gallery card", () => doubleClickOn(page, el("Model-Hub gallery card")));
    await session.step(30, "Then the tutorial step \"Run Model Hub\" should be done", () => stepDone(page, "Run Model Hub"));
    await session.step(31, "Given the tutorial step \"Run the \\\"Ball flight\\\" model\" should not be done yet", () => stepNotDone(page, "Run the \"Ball flight\" model"));
    await session.step(32, "When user double-clicks on \"Ball Flight Simulation\" link", () => doubleClickOn(page, el("\"Ball Flight Simulation\" link")));
    await session.step(33, "Then the tutorial step \"Run the \\\"Ball flight\\\" model\" should be done", () => stepDone(page, "Run the \"Ball flight\" model"));
    await session.step(34, "Given the tutorial step \"Click \\\"OK\\\"\" should not be done yet", () => stepNotDone(page, "Click \"OK\""));
    await session.step(35, "When user goes through the tour to its end", () => walkTour(page));
    await session.step(36, "Then the tutorial step \"Click \\\"OK\\\"\" should be done", () => stepDone(page, "Click \"OK\""));
    await session.step(38, "When user clicks on \"Fit inputs\" icon", () => clickOn(page, el("\"Fit inputs\" icon")));
    await session.step(39, "Then the tutorial step \"Click \\\"Fit inputs\\\"\" should be done", () => stepDone(page, "Click \"Fit inputs\""));
    await session.step(40, "Given the tutorial step \"Toggle the \\\"Velocity\\\" parameter\" should not be done yet", () => stepNotDone(page, "Toggle the \"Velocity\" parameter"));
    await session.step(41, "When user switches on \"Velocity\" input", () => switchOn(page, el("\"Velocity\" input")));
    await session.step(42, "Then the tutorial step \"Toggle the \\\"Velocity\\\" parameter\" should be done", () => stepDone(page, "Toggle the \"Velocity\" parameter"));
    await session.step(43, "When user switches on \"Angle\" input", () => switchOn(page, el("\"Angle\" input")));
    await session.step(44, "Then the tutorial step \"Toggle the \\\"Angle\\\" parameter\" should be done", () => stepDone(page, "Toggle the \"Angle\" parameter"));
    await session.step(45, "When user enters \"10\" into \"Max distance\" input", () => enterInto(page, "10", el("\"Max distance\" input")));
    await session.step(46, "Then the tutorial step \"Set \\\"Max distance\\\" to 10\" should be done", () => stepDone(page, "Set \"Max distance\" to 10"));
    await session.step(47, "When user clicks on \"Run\" icon", () => clickOn(page, el("\"Run\" icon")));
    await session.step(48, "Then the tutorial step \"Click \\\"Run\\\"\" should be done", () => stepDone(page, "Click \"Run\""));
    await session.step(49, "And bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
    await session.step(50, "Given the tutorial step \"Explore results\" should not be done yet", () => stepNotDone(page, "Explore results"));
    await session.step(51, "When user goes through the tour to its end", () => walkTour(page));
    await session.step(52, "Then the tutorial step \"Explore results\" should be done", () => stepDone(page, "Explore results"));
    await session.step(54, "Given the tutorial step \"Disable \\\"Max distance\\\"\" should not be done yet", () => stepNotDone(page, "Disable \"Max distance\""));
    await session.step(55, "When user switches off \"Max distance\" input", () => switchOff(page, el("\"Max distance\" input")));
    await session.step(56, "Then the tutorial step \"Disable \\\"Max distance\\\"\" should be done", () => stepDone(page, "Disable \"Max distance\""));
    await session.step(57, "When user switches on \"Trajectory\" input", () => switchOn(page, el("\"Trajectory\" input")));
    await session.step(58, "Then the tutorial step \"Toggle \\\"Trajectory\\\"\" should be done", () => stepDone(page, "Toggle \"Trajectory\""));
    await session.step(59, "When user selects \"Ball trajectory\" in \"Trajectory\" input", () => selectIn(page, "Ball trajectory", el("\"Trajectory\" input")));
    await session.step(60, "Then the tutorial step \"Set \\\"Trajectory\\\" to \\\"Ball trajectory\\\"\" should be done", () => stepDone(page, "Set \"Trajectory\" to \"Ball trajectory\""));
    await session.step(61, "When user clicks on \"Run\" icon", () => clickOn(page, el("\"Run\" icon")));
    await session.step(62, "Then the tutorial step \"Click \\\"Run\\\"\" should be done 2 times", () => stepDoneTimes(page, "Click \"Run\"", 2));
    await session.step(63, "Given the tutorial step \"Explore the fitted trajectory\" should not be done yet", () => stepNotDone(page, "Explore the fitted trajectory"));
    await session.step(64, "When user goes through the tour to its end", () => walkTour(page));
    await session.step(65, "Then the tutorial step \"Explore the fitted trajectory\" should be done", () => stepDone(page, "Explore the fitted trajectory"));
    await session.step(67, "And the \"Parameter optimization\" tutorial should be completed", () => tutorialCompleted(page, "Parameter optimization"));
    await session.step(68, "And the tutorial should have listed 15 steps", () => tutorialStepsListed(page, 15));
    await session.step(69, "And the tutorial progress should be 15 of 15", () => tutorialProgress(page, 15, 15));
    await session.step(70, "And no hint should be shown", () => noHintShown(page));
    await session.step(71, "And no errors should have been logged", () => noErrors(page));
  });
});
