/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/app/tutorials-app.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [apps.tutorials]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {cardDone, closeTutorial, startTutorial, stepNotDone, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialsClosed, tutorialsOpen} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, pressKeyIn, shouldBe, shouldContainText, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, dragSelectionOverArea, dragZoomOverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {doubleClickEmptySpace, pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Tutorials application", () => {
  const session = feature(test, "features/app/tutorials-app.feature", import.meta.url);
  test("Browse > Apps > Tutorials lists every track and tutorial", {tag: ["@tutorials", "@serial", "@realizes:apps.tutorials"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(17, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(20, "Given the Tutorials app is closed", () => tutorialsClosed(page));
    await session.step(21, "When user clicks on Apps tree node inside browse tree", () => clickOn(page, el("Apps tree node inside browse tree")));
    await session.step(22, "And user double-clicks on Tutorials gallery card", () => doubleClickOn(page, el("Tutorials gallery card")));
    await session.step(23, "Then Tutorials panel should be visible", () => shouldBe(page, el("Tutorials panel"), "visible"));
    await session.step(24, "And there should be 7 visible tutorial tracks", () => visibleCount(page, 7, el("tutorial tracks")));
    await session.step(25, "And there should be 20 visible tutorial cards", () => visibleCount(page, 20, el("tutorial cards")));
    await session.step(26, "And \"Exploratory data analysis\" tutorial track should be visible", () => shouldBe(page, el("\"Exploratory data analysis\" tutorial track"), "visible"));
    await session.step(27, "And \"Scatter Plot\" tutorial card should be visible", () => shouldBe(page, el("\"Scatter Plot\" tutorial card"), "visible"));
    await session.step(28, "And no errors should have been logged", () => noErrors(page));
  });
  test("A finished tutorial offers the next one of its track, which starts", {tag: ["@tutorials", "@serial", "@realizes:apps.tutorials"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(17, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(31, "Given the \"Scatter Plot\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Scatter Plot"));
    await session.step(32, "And the \"Embedded Viewers\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Embedded Viewers"));
    await session.step(33, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(34, "When user starts the \"Scatter Plot\" tutorial", () => startTutorial(page, "Scatter Plot"));
    await session.step(35, "And user clicks on scatter-plot icon in toolbox", () => clickOn(page, el("scatter-plot icon in toolbox")));
    await session.step(36, "And user picks \"HEIGHT\" in the \"x\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "HEIGHT", "x", el("scatter plot viewer")));
    await session.step(37, "And user picks \"WEIGHT\" in the \"y\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "WEIGHT", "y", el("scatter plot viewer")));
    await session.step(38, "And user picks \"AGE\" in the \"size\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "AGE", "size", el("scatter plot viewer")));
    await session.step(39, "And user picks \"SEX\" in the \"color\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "SEX", "color", el("scatter plot viewer")));
    await session.step(40, "And user drags a zoom box over the \"view\" area of scatter plot viewer", () => dragZoomOverArea(page, "view", el("scatter plot viewer")));
    await session.step(41, "And user double-clicks on empty plot space of scatter plot viewer", () => doubleClickEmptySpace(page, el("scatter plot viewer")));
    await session.step(42, "And user clicks on the \"marker of row 11\" area of scatter plot viewer", () => clickArea(page, "marker of row 11", el("scatter plot viewer")));
    await session.step(43, "And user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
    await session.step(44, "And user presses Escape in scatter plot viewer", () => pressKeyIn(page, "Escape", el("scatter plot viewer")));
    await session.step(45, "Then the \"Scatter Plot\" tutorial should be completed", () => tutorialCompleted(page, "Scatter Plot"));
    await session.step(46, "And Tutorials panel should contain text \"Next \\\"Embedded Viewers\\\"\"", () => shouldContainText(page, el("Tutorials panel"), "Next \"Embedded Viewers\""));
    await session.step(48, "When user clicks on Start button in Tutorials panel", () => clickOn(page, el("Start button in Tutorials panel")));
    await session.step(49, "Then tutorial title should contain text \"Embedded Viewers\"", () => shouldContainText(page, el("tutorial title"), "Embedded Viewers"));
    await session.step(50, "And the tutorial progress should be 1 of 10", () => tutorialProgress(page, 1, 10));
    await session.step(51, "And the tutorial step \"Open scatter plot\" should not be done yet", () => stepNotDone(page, "Open scatter plot"));
    await session.step(52, "When user closes the tutorial", () => closeTutorial(page));
    await session.step(53, "Then the \"Scatter Plot\" tutorial card should show it is done", () => cardDone(page, "Scatter Plot"));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
  });
});
