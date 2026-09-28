/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/compute/differential-equations.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.differential-equations]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {startTutorial, stepDone, stepNotDone, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, insertLineAfter, replaceLine, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, noHintShown, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Differential Equations tutorial", () => {
  const session = feature(test, "features/compute/differential-equations.feature", import.meta.url);
  test("A learner completes the Differential equations tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.differential-equations"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(24, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(25, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(26, "And the \"Differential equations\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Differential equations"));
    await session.step(27, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(30, "When user starts the \"Differential equations\" tutorial", () => startTutorial(page, "Differential equations"));
    await session.step(31, "Then the tutorial progress should be 1 of 14", () => tutorialProgress(page, 1, 14));
    await session.step(32, "When user clicks on Apps tree node inside browse tree", () => clickOn(page, el("Apps tree node inside browse tree")));
    await session.step(33, "Then the tutorial step \"Open Apps\" should be done", () => stepDone(page, "Open Apps"));
    await session.step(35, "Given the tutorial step \"Run Diff Studio\" should not be done yet", () => stepNotDone(page, "Run Diff Studio"));
    await session.step(36, "When user double-clicks on Diff-Studio gallery card", () => doubleClickOn(page, el("Diff-Studio gallery card")));
    await session.step(37, "Then the tutorial step \"Run Diff Studio\" should be done", () => stepDone(page, "Run Diff Studio"));
    await session.step(38, "Given the tutorial step \"Run the Lotka-Volterra model\" should not be done yet", () => stepNotDone(page, "Run the Lotka-Volterra model"));
    await session.step(39, "When user double-clicks on \"Lotka-Volterra\" model card", () => doubleClickOn(page, el("\"Lotka-Volterra\" model card")));
    await session.step(40, "Then the tutorial step \"Run the Lotka-Volterra model\" should be done", () => stepDone(page, "Run the Lotka-Volterra model"));
    await session.step(42, "Given the tutorial step \"Explore the interface\" should not be done yet", () => stepNotDone(page, "Explore the interface"));
    await session.step(43, "When user clicks on \"next\" tour button", () => clickOn(page, el("\"next\" tour button")));
    await session.step(44, "And user clicks on \"next\" tour button", () => clickOn(page, el("\"next\" tour button")));
    await session.step(45, "And user clicks on \"next\" tour button", () => clickOn(page, el("\"next\" tour button")));
    await session.step(46, "And user clicks on \"next\" tour button", () => clickOn(page, el("\"next\" tour button")));
    await session.step(47, "And user clicks on \"done\" tour button", () => clickOn(page, el("\"done\" tour button")));
    await session.step(48, "Then the tutorial step \"Explore the interface\" should be done", () => stepDone(page, "Explore the interface"));
    await session.step(50, "Given the tutorial step \"Open equations editor\" should not be done yet", () => stepNotDone(page, "Open equations editor"));
    await session.step(51, "When user clicks on Edit ribbon item", () => clickOn(page, el("Edit ribbon item")));
    await session.step(52, "Then the tutorial step \"Open equations editor\" should be done", () => stepDone(page, "Open equations editor"));
    await session.step(53, "And code editor should be visible", () => shouldBe(page, el("code editor"), "visible"));
    await session.step(54, "Given the tutorial step \"Explore editor\" should not be done yet", () => stepNotDone(page, "Explore editor"));
    await session.step(55, "When user clicks on \"next\" tour button", () => clickOn(page, el("\"next\" tour button")));
    await session.step(56, "And user clicks on \"next\" tour button", () => clickOn(page, el("\"next\" tour button")));
    await session.step(57, "And user clicks on \"next\" tour button", () => clickOn(page, el("\"next\" tour button")));
    await session.step(58, "And user clicks on \"next\" tour button", () => clickOn(page, el("\"next\" tour button")));
    await session.step(59, "And user clicks on \"next\" tour button", () => clickOn(page, el("\"next\" tour button")));
    await session.step(60, "And user clicks on \"done\" tour button", () => clickOn(page, el("\"done\" tour button")));
    await session.step(61, "Then the tutorial step \"Explore editor\" should be done", () => stepDone(page, "Explore editor"));
    await session.step(63, "Given the tutorial step \"Complete the predator equation\" should not be done yet", () => stepNotDone(page, "Complete the predator equation"));
    await session.step(64, "When user replaces the line starting with \"dy/dt\" in code editor with \"dy/dt = -gamma * y + delta * x * y - eta * y * y\"", () => replaceLine(page, "dy/dt", el("code editor"), "dy/dt = -gamma * y + delta * x * y - eta * y * y"));
    await session.step(65, "Then the tutorial step \"Complete the predator equation\" should be done", () => stepDone(page, "Complete the predator equation"));
    await session.step(66, "When user puts \"eta = 0.01 {min: 0; max: 0.1; category: Parameters} [Crowding effect]\" on a new line after the line starting with \"#parameters\" in code editor", () => insertLineAfter(page, "eta = 0.01 {min: 0; max: 0.1; category: Parameters} [Crowding effect]", "#parameters", el("code editor")));
    await session.step(67, "Then the tutorial step \"Add the eta parameter\" should be done", () => stepDone(page, "Add the eta parameter"));
    await session.step(68, "When user clicks on Refresh ribbon item", () => clickOn(page, el("Refresh ribbon item")));
    await session.step(69, "Then the tutorial step \"Apply changes\" should be done", () => stepDone(page, "Apply changes"));
    await session.step(70, "Given the tutorial step \"Check the updates\" should not be done yet", () => stepNotDone(page, "Check the updates"));
    await session.step(71, "When user clicks on \"done\" tour button", () => clickOn(page, el("\"done\" tour button")));
    await session.step(72, "Then the tutorial step \"Check the updates\" should be done", () => stepDone(page, "Check the updates"));
    await session.step(74, "When user clicks on Edit ribbon item", () => clickOn(page, el("Edit ribbon item")));
    await session.step(75, "Then the tutorial step \"Close equations editor\" should be done", () => stepDone(page, "Close equations editor"));
    await session.step(76, "And code editor should be absent", () => shouldBe(page, el("code editor"), "absent"));
    await session.step(77, "When user enters \"2\" into \"Prey\" input", () => enterInto(page, "2", el("\"Prey\" input")));
    await session.step(78, "Then the tutorial step \"Set \\\"Prey\\\" to 2\" should be done", () => stepDone(page, "Set \"Prey\" to 2"));
    await session.step(79, "When user enters \"0.1\" into \"Delta\" input", () => enterInto(page, "0.1", el("\"Delta\" input")));
    await session.step(80, "Then the tutorial step \"Set \\\"Delta\\\" to 0.1\" should be done", () => stepDone(page, "Set \"Delta\" to 0.1"));
    await session.step(81, "When user enters \"150\" into \"Finish\" input", () => enterInto(page, "150", el("\"Finish\" input")));
    await session.step(82, "Then the tutorial step \"Set \\\"Finish\\\" to 150\" should be done", () => stepDone(page, "Set \"Finish\" to 150"));
    await session.step(84, "And the \"Differential equations\" tutorial should be completed", () => tutorialCompleted(page, "Differential equations"));
    await session.step(85, "And the tutorial should have listed 14 steps", () => tutorialStepsListed(page, 14));
    await session.step(86, "And the tutorial progress should be 14 of 14", () => tutorialProgress(page, 14, 14));
    await session.step(87, "And no hint should be shown", () => noHintShown(page));
    await session.step(88, "And no errors should have been logged", () => noErrors(page));
  });
});
