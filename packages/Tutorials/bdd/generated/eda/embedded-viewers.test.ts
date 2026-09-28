/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/eda/embedded-viewers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.embedded-viewers]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {startTutorial, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {elementHinted, noHintShown, toolboxPaneHidden, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noErrors, pickFromAreaContextMenu, pickFromContextMenu, readingNotAsRemembered, readingReads, rememberReading, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {innerPropertyShouldBe, pickInColumnSelector, pickInnerViewer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Embedded Viewers tutorial", () => {
  const session = feature(test, "features/eda/embedded-viewers.feature", import.meta.url);
  test("A learner completes the Embedded Viewers tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.embedded-viewers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(21, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(22, "And the \"recentViewerSettings\" user settings are put back at feature end", () => userSettingsPutBack(page, "recentViewerSettings"));
    await session.step(23, "And the \"Embedded Viewers\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Embedded Viewers"));
    await session.step(24, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(27, "When user starts the \"Embedded Viewers\" tutorial", () => startTutorial(page, "Embedded Viewers"));
    await session.step(28, "Then the tutorial progress should be 1 of 10", () => tutorialProgress(page, 1, 10));
    await session.step(29, "When user clicks on scatter-plot icon in toolbox", () => clickOn(page, el("scatter-plot icon in toolbox")));
    await session.step(30, "Then \"Open scatter plot\" tutorial step should be checked", () => shouldBe(page, el("\"Open scatter plot\" tutorial step"), "checked"));
    await session.step(31, "When user clicks on histogram icon in toolbox", () => clickOn(page, el("histogram icon in toolbox")));
    await session.step(32, "Then \"Open histogram\" tutorial step should be checked", () => shouldBe(page, el("\"Open histogram\" tutorial step"), "checked"));
    await session.step(33, "When user clicks on pie-chart icon in toolbox", () => clickOn(page, el("pie-chart icon in toolbox")));
    await session.step(34, "Then \"Open pie chart\" tutorial step should be checked", () => shouldBe(page, el("\"Open pie chart\" tutorial step"), "checked"));
    await session.step(35, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
    await session.step(36, "And the open tableview should have 1 histogram viewer", () => viewerCount(page, 1, "histogram"));
    await session.step(37, "And the open tableview should have 1 pie chart viewer", () => viewerCount(page, 1, "pie chart"));
    await session.step(39, "When user picks \"Tooltip > Use as Group Tooltip\" from the context menu of the \"view\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Tooltip > Use as Group Tooltip", "view", el("scatter plot viewer")));
    await session.step(40, "Then \"Set the scatter plot as a tooltip viewer\" tutorial step should be checked", () => shouldBe(page, el("\"Set the scatter plot as a tooltip viewer\" tutorial step"), "checked"));
    await session.step(41, "When user hovers over the \"bin 3\" area of histogram viewer", () => hoverArea(page, "bin 3", el("histogram viewer")));
    await session.step(42, "Then \"Hover over the histogram bins or pie chart segments\" tutorial step should be checked", () => shouldBe(page, el("\"Hover over the histogram bins or pie chart segments\" tutorial step"), "checked"));
    await session.step(43, "And scatter plot viewer in tooltip should be visible", () => shouldBe(page, el("scatter plot viewer in tooltip"), "visible"));
    await session.step(45, "When user picks \"Tooltip > Remove Group Tooltip\" from the context menu of the \"view\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Tooltip > Remove Group Tooltip", "view", el("scatter plot viewer")));
    await session.step(46, "Then \"Reset the tooltip\" tutorial step should be checked", () => shouldBe(page, el("\"Reset the tooltip\" tutorial step"), "checked"));
    await session.step(47, "When user hovers over the \"bin 4\" area of histogram viewer", () => hoverArea(page, "bin 4", el("histogram viewer")));
    await session.step(48, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
    await session.step(49, "And scatter plot viewer in tooltip should be absent", () => shouldBe(page, el("scatter plot viewer in tooltip"), "absent"));
    await session.step(51, "When user picks \"General > Use in Trellis\" from the context menu of pie chart viewer", () => pickFromContextMenu(page, "General > Use in Trellis", el("pie chart viewer")));
    await session.step(52, "Then \"Open a Trellis plot from the pie chart's context menu\" tutorial step should be checked", () => shouldBe(page, el("\"Open a Trellis plot from the pie chart's context menu\" tutorial step"), "checked"));
    await session.step(53, "And the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
    await session.step(54, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Pie chart\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Pie chart"));
    await session.step(56, "Then viewer selector in trellis plot viewer should be hinted", () => elementHinted(page, el("viewer selector in trellis plot viewer")));
    await session.step(57, "When user picks \"Scatter plot\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Scatter plot", el("trellis plot viewer")));
    await session.step(58, "Then \"Set a scatter plot as an inner viewer\" tutorial step should be checked", () => shouldBe(page, el("\"Set a scatter plot as an inner viewer\" tutorial step"), "checked"));
    await session.step(59, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Scatter plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Scatter plot"));
    await session.step(64, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(65, "When user presses F4", () => pressKey(page, "F4"));
    await session.step(66, "Then context panel should be hidden", () => shouldBe(page, el("context panel"), "hidden"));
    await session.step(67, "And color column selector in trellis plot viewer should be hinted", () => elementHinted(page, el("color column selector in trellis plot viewer")));
    await session.step(70, "When user remembers the \"cell signature RA | Critical\" reading of trellis plot viewer", () => rememberReading(page, "cell signature RA | Critical", el("trellis plot viewer")));
    await session.step(71, "And user picks \"AGE\" in the \"color\" column selector of trellis plot viewer", () => pickInColumnSelector(page, "AGE", "color", el("trellis plot viewer")));
    await session.step(72, "Then \"Set Color of the inner scatter plot to AGE\" tutorial step should be checked", () => shouldBe(page, el("\"Set Color of the inner scatter plot to AGE\" tutorial step"), "checked"));
    await session.step(73, "And \"colorColumnName\" inner property of trellis plot viewer should be \"AGE\"", () => innerPropertyShouldBe(page, "colorColumnName", el("trellis plot viewer"), "AGE"));
    await session.step(74, "And the \"cell signature RA | Critical\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature RA | Critical", el("trellis plot viewer")));
    await session.step(76, "And the \"Embedded Viewers\" tutorial should be completed", () => tutorialCompleted(page, "Embedded Viewers"));
    await session.step(77, "And the tutorial should have listed 9 steps", () => tutorialStepsListed(page, 9));
    await session.step(78, "And the tutorial progress should be 10 of 10", () => tutorialProgress(page, 10, 10));
    await session.step(79, "And no hint should be shown", () => noHintShown(page));
    await session.step(80, "And no errors should have been logged", () => noErrors(page));
  });
});
