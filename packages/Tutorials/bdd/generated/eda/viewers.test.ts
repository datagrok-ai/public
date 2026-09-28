/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/eda/viewers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.viewers]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {startTutorial, stepDone, stepDoneTimes, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, rememberVisibleCount, shouldBe, typeInto, visibleFewerThanRemembered} from '@datagrok-libraries/bdd/bindings/common/steps';
import {onlyOfSelected, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, contextPanelShows, elementHinted, noHintShown, userSettingsPutBack, viewHoldsViewers} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, dragSelectionOverArea, hoverArea, noErrors, pickFromAreaContextMenu, propertyShouldBe, readingDoesNotRead, readingNotAsRemembered, rememberReading, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {currentRowIsHovered, pickFromViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Viewers tutorial", () => {
  const session = feature(test, "features/eda/viewers.feature", import.meta.url);
  test("A learner completes the Viewers tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.viewers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(22, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(23, "And the \"recentViewerSettings\" user settings are put back at feature end", () => userSettingsPutBack(page, "recentViewerSettings"));
    await session.step(24, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(25, "And the \"Viewers\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Viewers"));
    await session.step(26, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(29, "When user starts the \"Viewers\" tutorial", () => startTutorial(page, "Viewers"));
    await session.step(30, "Then the tutorial progress should be 1 of 21", () => tutorialProgress(page, 1, 21));
    await session.step(31, "And \"Add viewer\" icon should be hinted", () => elementHinted(page, el("\"Add viewer\" icon")));
    await session.step(32, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(33, "Then viewer gallery should be visible", () => shouldBe(page, el("viewer gallery"), "visible"));
    await session.step(34, "And the tutorial step \"Click the Add viewer icon to open the gallery\" should be done", () => stepDone(page, "Click the Add viewer icon to open the gallery"));
    await session.step(36, "Then \"Charts\" viewer tag should be hinted", () => elementHinted(page, el("\"Charts\" viewer tag")));
    await session.step(37, "When user remembers the number of visible viewer cards", () => rememberVisibleCount(page, el("viewer cards")));
    await session.step(38, "And user clicks on \"Charts\" viewer tag", () => clickOn(page, el("\"Charts\" viewer tag")));
    await session.step(39, "Then the tutorial step \"Click the \\\"Charts\\\" tag to filter the viewers\" should be done", () => stepDone(page, "Click the \"Charts\" tag to filter the viewers"));
    await session.step(40, "And there should be fewer visible viewer cards than remembered", () => visibleFewerThanRemembered(page, el("viewer cards")));
    await session.step(41, "When user clicks on \"Radar\" viewer card", () => clickOn(page, el("\"Radar\" viewer card")));
    await session.step(42, "Then the tutorial step \"Select the Radar viewer\" should be done", () => stepDone(page, "Select the Radar viewer"));
    await session.step(43, "And the current view should hold at least 2 viewers", () => viewHoldsViewers(page, 2));
    await session.step(45, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(46, "Then the tutorial step \"Open the viewer gallery again\" should be done", () => stepDone(page, "Open the viewer gallery again"));
    await session.step(47, "When user types \"Sunburst\" into viewer gallery search", () => typeInto(page, "Sunburst", el("viewer gallery search")));
    await session.step(48, "Then the tutorial step \"Type \\\"Sunburst\\\" in the search box\" should be done", () => stepDone(page, "Type \"Sunburst\" in the search box"));
    await session.step(49, "When user clicks on \"Sunburst\" viewer card", () => clickOn(page, el("\"Sunburst\" viewer card")));
    await session.step(50, "Then the tutorial step \"Select the Sunburst viewer\" should be done", () => stepDone(page, "Select the Sunburst viewer"));
    await session.step(51, "And the open tableview should have 1 sunburst viewer", () => viewerCount(page, 1, "sunburst"));
    await session.step(53, "When user clicks on scatter-plot icon in toolbox", () => clickOn(page, el("scatter-plot icon in toolbox")));
    await session.step(54, "Then the tutorial step \"Open scatter plot\" should be done", () => stepDone(page, "Open scatter plot"));
    await session.step(55, "When user clicks on histogram icon in toolbox", () => clickOn(page, el("histogram icon in toolbox")));
    await session.step(56, "Then the tutorial step \"Open histogram\" should be done", () => stepDone(page, "Open histogram"));
    await session.step(57, "When user clicks on pie-chart icon in toolbox", () => clickOn(page, el("pie-chart icon in toolbox")));
    await session.step(58, "Then the tutorial step \"Open pie chart\" should be done", () => stepDone(page, "Open pie chart"));
    await session.step(60, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
    await session.step(61, "Then the tutorial step \"Hover over the histogram bins or scatter plot points\" should be done", () => stepDone(page, "Hover over the histogram bins or scatter plot points"));
    await session.step(62, "And the \"hovered row\" reading of scatter plot viewer should not be \"0\"", () => readingDoesNotRead(page, "hovered row", el("scatter plot viewer"), "0"));
    await session.step(64, "When user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
    await session.step(65, "Then the tutorial step \"Select points on the scatter plot\" should be done", () => stepDone(page, "Select points on the scatter plot"));
    await session.step(66, "And some rows should be selected", () => someSelected(page));
    await session.step(67, "When user remembers the \"rows selected\" reading of scatter plot viewer", () => rememberReading(page, "rows selected", el("scatter plot viewer")));
    await session.step(68, "And user clicks on the \"bin 3\" area of histogram viewer", () => clickArea(page, "bin 3", el("histogram viewer")));
    await session.step(69, "Then the tutorial step \"Select one of the bins on the histogram\" should be done", () => stepDone(page, "Select one of the bins on the histogram"));
    await session.step(70, "And the \"rows selected\" reading of scatter plot viewer should not be as remembered", () => readingNotAsRemembered(page, "rows selected", el("scatter plot viewer")));
    await session.step(71, "When user remembers the \"rows selected\" reading of scatter plot viewer", () => rememberReading(page, "rows selected", el("scatter plot viewer")));
    await session.step(73, "And user clicks on the \"segment M\" area of sunburst viewer", () => clickArea(page, "segment M", el("sunburst viewer")));
    await session.step(74, "Then the tutorial step \"Click a Sunburst segment to select its rows\" should be done", () => stepDone(page, "Click a Sunburst segment to select its rows"));
    await session.step(75, "And the \"rows selected\" reading of scatter plot viewer should not be as remembered", () => readingNotAsRemembered(page, "rows selected", el("scatter plot viewer")));
    await session.step(76, "And only rows where \"SEX\" is \"M\" should be selected", () => onlyOfSelected(page, "SEX", "M"));
    await session.step(78, "When user clicks on the \"marker of row 11\" area of scatter plot viewer", () => clickArea(page, "marker of row 11", el("scatter plot viewer")));
    await session.step(79, "Then the tutorial step \"Click on a point to set the current record\" should be done", () => stepDone(page, "Click on a point to set the current record"));
    await session.step(80, "And the \"hovered row\" reading of scatter plot viewer should be the current row", () => currentRowIsHovered(page, "hovered row", el("scatter plot viewer")));
    await session.step(82, "When user picks \"Properties...\" from the viewer menu of scatter plot viewer", () => pickFromViewerMenu(page, "Properties...", el("scatter plot viewer")));
    await session.step(83, "Then the tutorial step \"Open the scatter plot's properties\" should be done", () => stepDone(page, "Open the scatter plot's properties"));
    await session.step(84, "And the context panel should show \"Scatter plot\"", () => contextPanelShows(page, "Scatter plot"));
    await session.step(85, "When user sets \"markerDefaultSize\" property of scatter plot viewer to \"12\"", () => setProperty(page, "markerDefaultSize", el("scatter plot viewer"), "12"));
    await session.step(86, "Then the tutorial step \"Change a few visual properties, e.g., the background color or marker size\" should be done", () => stepDone(page, "Change a few visual properties, e.g., the background color or marker size"));
    await session.step(88, "When user picks \"General > Clone\" from the context menu of the \"view\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "General > Clone", "view", el("scatter plot viewer")));
    await session.step(89, "Then the tutorial step \"Clone the scatter plot\" should be done", () => stepDone(page, "Clone the scatter plot"));
    await session.step(90, "And the open tableview should have 2 scatter plot viewers", () => viewerCount(page, 2, "scatter plot"));
    await session.step(91, "When user picks \"Pick Up / Apply > Pick Up\" from the context menu of the \"view\" area of first scatter plot viewer", () => pickFromAreaContextMenu(page, "Pick Up / Apply > Pick Up", "view", el("first scatter plot viewer")));
    await session.step(92, "Then the tutorial step \"Pick up the scatter plot's style\" should be done", () => stepDone(page, "Pick up the scatter plot's style"));
    await session.step(93, "When user clicks on scatter-plot icon in toolbox", () => clickOn(page, el("scatter-plot icon in toolbox")));
    await session.step(94, "Then the tutorial step \"Open scatter plot\" should be done 2 times", () => stepDoneTimes(page, "Open scatter plot", 2));
    await session.step(95, "And the open tableview should have 3 scatter plot viewers", () => viewerCount(page, 3, "scatter plot"));
    await session.step(96, "When user picks \"Pick Up / Apply > Apply\" from the context menu of the \"view\" area of third scatter plot viewer", () => pickFromAreaContextMenu(page, "Pick Up / Apply > Apply", "view", el("third scatter plot viewer")));
    await session.step(97, "Then the tutorial step \"Apply the style to the new viewer\" should be done", () => stepDone(page, "Apply the style to the new viewer"));
    await session.step(98, "And \"markerDefaultSize\" property of third scatter plot viewer should be \"12\"", () => propertyShouldBe(page, "markerDefaultSize", el("third scatter plot viewer"), "12"));
    await session.step(100, "And the \"Viewers\" tutorial should be completed", () => tutorialCompleted(page, "Viewers"));
    await session.step(101, "And the tutorial should have listed 20 steps", () => tutorialStepsListed(page, 20));
    await session.step(102, "And the tutorial progress should be 21 of 21", () => tutorialProgress(page, 21, 21));
    await session.step(103, "And no hint should be shown", () => noHintShown(page));
    await session.step(104, "And no errors should have been logged", () => noErrors(page));
  });
});
