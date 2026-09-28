/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/eda/scatter-plot.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.scatter-plot]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {startTutorial, stepDone, stepNotDone, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKeyIn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {noneSelected, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {elementHinted, noHintShown, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, dragSelectionOverArea, dragZoomOverArea, noErrors, propertyShouldBe, readingAsRemembered, readingLowerThanRemembered, rememberReading, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {currentRowIsHovered, doubleClickEmptySpace, pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Scatter Plot tutorial", () => {
  const session = feature(test, "features/eda/scatter-plot.feature", import.meta.url);
  test("A learner completes the Scatter Plot tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(17, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(18, "And the \"recentViewerSettings\" user settings are put back at feature end", () => userSettingsPutBack(page, "recentViewerSettings"));
    await session.step(19, "And the \"Scatter Plot\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Scatter Plot"));
    await session.step(20, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(23, "When user starts the \"Scatter Plot\" tutorial", () => startTutorial(page, "Scatter Plot"));
    await session.step(24, "Then the tutorial progress should be 1 of 11", () => tutorialProgress(page, 1, 11));
    await session.step(25, "And the tutorial step \"Open scatter plot\" should not be done yet", () => stepNotDone(page, "Open scatter plot"));
    await session.step(26, "And scatter-plot icon in toolbox should be hinted", () => elementHinted(page, el("scatter-plot icon in toolbox")));
    await session.step(27, "When user clicks on scatter-plot icon in toolbox", () => clickOn(page, el("scatter-plot icon in toolbox")));
    await session.step(28, "Then the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
    await session.step(29, "And the tutorial step \"Open scatter plot\" should be done", () => stepDone(page, "Open scatter plot"));
    await session.step(30, "And the tutorial progress should be 2 of 11", () => tutorialProgress(page, 2, 11));
    await session.step(32, "When user picks \"HEIGHT\" in the \"x\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "HEIGHT", "x", el("scatter plot viewer")));
    await session.step(33, "Then the tutorial step \"Set X to HEIGHT\" should be done", () => stepDone(page, "Set X to HEIGHT"));
    await session.step(34, "And \"xColumnName\" property of scatter plot viewer should be \"HEIGHT\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer"), "HEIGHT"));
    await session.step(35, "When user picks \"WEIGHT\" in the \"y\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "WEIGHT", "y", el("scatter plot viewer")));
    await session.step(36, "Then the tutorial step \"Set Y to WEIGHT\" should be done", () => stepDone(page, "Set Y to WEIGHT"));
    await session.step(37, "And \"yColumnName\" property of scatter plot viewer should be \"WEIGHT\"", () => propertyShouldBe(page, "yColumnName", el("scatter plot viewer"), "WEIGHT"));
    await session.step(38, "When user picks \"AGE\" in the \"size\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "AGE", "size", el("scatter plot viewer")));
    await session.step(39, "Then the tutorial step \"Set Size to AGE\" should be done", () => stepDone(page, "Set Size to AGE"));
    await session.step(40, "And \"sizeColumnName\" property of scatter plot viewer should be \"AGE\"", () => propertyShouldBe(page, "sizeColumnName", el("scatter plot viewer"), "AGE"));
    await session.step(41, "When user picks \"SEX\" in the \"color\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "SEX", "color", el("scatter plot viewer")));
    await session.step(42, "Then the tutorial step \"Set Color to SEX\" should be done", () => stepDone(page, "Set Color to SEX"));
    await session.step(43, "And \"colorColumnName\" property of scatter plot viewer should be \"SEX\"", () => propertyShouldBe(page, "colorColumnName", el("scatter plot viewer"), "SEX"));
    await session.step(45, "When user remembers the \"x axis span\" reading of scatter plot viewer", () => rememberReading(page, "x axis span", el("scatter plot viewer")));
    await session.step(46, "And user drags a zoom box over the \"view\" area of scatter plot viewer", () => dragZoomOverArea(page, "view", el("scatter plot viewer")));
    await session.step(47, "Then the tutorial step \"Zoom in\" should be done", () => stepDone(page, "Zoom in"));
    await session.step(48, "And the \"x axis span\" reading of scatter plot viewer should be lower than remembered", () => readingLowerThanRemembered(page, "x axis span", el("scatter plot viewer")));
    await session.step(49, "When user double-clicks on empty plot space of scatter plot viewer", () => doubleClickEmptySpace(page, el("scatter plot viewer")));
    await session.step(50, "Then the tutorial step \"Double-click to unzoom\" should be done", () => stepDone(page, "Double-click to unzoom"));
    await session.step(51, "And the \"x axis span\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "x axis span", el("scatter plot viewer")));
    await session.step(55, "When user clicks on the \"marker of row 11\" area of scatter plot viewer", () => clickArea(page, "marker of row 11", el("scatter plot viewer")));
    await session.step(56, "Then the tutorial step \"Click on a point\" should be done", () => stepDone(page, "Click on a point"));
    await session.step(57, "And the \"hovered row\" reading of scatter plot viewer should be the current row", () => currentRowIsHovered(page, "hovered row", el("scatter plot viewer")));
    await session.step(59, "When user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
    await session.step(60, "Then the tutorial step \"Select points\" should be done", () => stepDone(page, "Select points"));
    await session.step(61, "And some rows should be selected", () => someSelected(page));
    await session.step(62, "And the tutorial step \"Deselect points\" should not be done yet", () => stepNotDone(page, "Deselect points"));
    await session.step(63, "When user presses Escape in scatter plot viewer", () => pressKeyIn(page, "Escape", el("scatter plot viewer")));
    await session.step(64, "Then the tutorial step \"Deselect points\" should be done", () => stepDone(page, "Deselect points"));
    await session.step(65, "And no rows should be selected", () => noneSelected(page));
    await session.step(67, "And the \"Scatter Plot\" tutorial should be completed", () => tutorialCompleted(page, "Scatter Plot"));
    await session.step(68, "And the tutorial should have listed 10 steps", () => tutorialStepsListed(page, 10));
    await session.step(69, "And the tutorial progress should be 11 of 11", () => tutorialProgress(page, 11, 11));
    await session.step(70, "And no hint should be shown", () => noHintShown(page));
    await session.step(71, "And no errors should have been logged", () => noErrors(page));
  });
});
