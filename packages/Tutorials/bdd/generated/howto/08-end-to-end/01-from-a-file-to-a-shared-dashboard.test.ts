/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/08-end-to-end/01-from-a-file-to-a-shared-dashboard.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, collapse, doubleClickOn, enterInto, expand, hoverOver, shouldBe, shouldContainText, typeInto, uncheck, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterIsExactlyCategory, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, contextPanelClosed, contextPanelShows, dialogCloses, noProjectOnServer, projectsOnServer, sharingUserNoNotifications, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, infoBalloonText, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · End to end", () => {
  const session = feature(test, "features/howto/08-end-to-end/01-from-a-file-to-a-shared-dashboard.feature", import.meta.url);
  test("From a file to a shared dashboard", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
    await session.step(21, "And no project named \"Rheumatoid arthritis overview\" is on the server", () => noProjectOnServer(page, "Rheumatoid arthritis overview"));
    await session.step(22, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(23, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "When the video says \"The task: a chart overview of one patient group, shared with a colleague\"", () => videoSays(page, "The task: a chart overview of one patient group, shared with a colleague"));
    await session.step(27, "And user expands Files tree node inside browse tree", () => expand(page, el("Files tree node inside browse tree")), undefined, "Start from the file: Browse > Files > Demo holds sample data");
    await session.step(28, "And user expands Files---Demo tree node inside browse tree", () => expand(page, el("Files---Demo tree node inside browse tree")));
    await session.step(30, "And user double-clicks on Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")), undefined, "Double-click demog.csv to open it");
    await session.step(31, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(32, "And the table should have 5850 rows", () => rowCount(page, 5850));
    await session.step(34, "When user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")), undefined, "First, look: switch on the Context Panel and click a column");
    await session.step(35, "And user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")));
    await session.step(36, "Then the context panel should show \"AGE\"", () => contextPanelShows(page, "AGE"));
    await session.step(38, "When user expands Stats pane in context panel", () => expand(page, el("Stats pane in context panel")), undefined, "Stats sums the column up");
    await session.step(39, "Then Stats pane in context panel should contain text \"avg\"", () => shouldContainText(page, el("Stats pane in context panel"), "avg"));
    await session.step(40, "When user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")));
    await session.step(41, "Then toolbox should be visible", () => shouldBe(page, el("toolbox"), "visible"));
    await session.step(43, "When user clicks on \"bar chart\" icon in toolbox", () => clickOn(page, el("\"bar chart\" icon in toolbox")), undefined, "Chart it from the toolbox: a bar chart counts the patients per group");
    await session.step(44, "Then bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
    await session.step(45, "When user picks \"DIS_POP\" in the \"split\" column selector of bar chart viewer", () => pickInColumnSelector(page, "DIS_POP", "split", el("bar chart viewer")));
    await session.step(47, "And user clicks on \"scatter plot\" icon in toolbox", () => clickOn(page, el("\"scatter plot\" icon in toolbox")), undefined, "A scatter plot: weight against height, colored by sex");
    await session.step(48, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(49, "When user picks \"HEIGHT\" in the \"x\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "HEIGHT", "x", el("scatter plot viewer")));
    await session.step(50, "And user picks \"WEIGHT\" in the \"y\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "WEIGHT", "y", el("scatter plot viewer")));
    await session.step(51, "And user picks \"SEX\" in the \"color\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "SEX", "color", el("scatter plot viewer")));
    await session.step(53, "And user clicks on \"histogram\" icon in toolbox", () => clickOn(page, el("\"histogram\" icon in toolbox")), undefined, "A histogram shows the ages");
    await session.step(54, "Then histogram viewer should be visible", () => shouldBe(page, el("histogram viewer"), "visible"));
    await session.step(55, "When user picks \"AGE\" in the \"value\" column selector of histogram viewer", () => pickInColumnSelector(page, "AGE", "value", el("histogram viewer")));
    await session.step(57, "And user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")), undefined, "Narrow it to the group that matters: the filter panel");
    await session.step(58, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
    await session.step(59, "When user clicks on the \"category RA of DIS_POP\" area of filter panel", () => clickArea(page, "category RA of DIS_POP", el("filter panel")));
    await session.step(60, "Then the filter should pass exactly the rows where \"DIS_POP\" is \"RA\"", () => filterIsExactlyCategory(page, "DIS_POP", "RA"));
    await session.step(61, "When the video says \"Every chart now shows only the RA patients\"", () => videoSays(page, "Every chart now shows only the RA patients"));
    await session.step(62, "And user hovers over scatter plot viewer", () => hoverOver(page, el("scatter plot viewer")));
    await session.step(64, "And user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")), undefined, "SAVE keeps the table, the charts and the filter together");
    await session.step(65, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(66, "When user enters \"Rheumatoid arthritis overview\" into Name text input in \"Save project\" dialog", () => enterInto(page, "Rheumatoid arthritis overview", el("Name text input in \"Save project\" dialog")));
    await session.step(67, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(68, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(69, "And 1 project named \"Rheumatoid arthritis overview\" should be on the server", () => projectsOnServer(page, 1, "Rheumatoid arthritis overview"));
    await session.step(70, "And \"Share Rheumatoid arthritis overview\" dialog should be visible", () => shouldBe(page, el("\"Share Rheumatoid arthritis overview\" dialog"), "visible"));
    await session.step(72, "When user types \"bddsecond\" into \"User, group, or email\" input in \"Share Rheumatoid arthritis overview\" dialog", () => typeInto(page, "bddsecond", el("\"User, group, or email\" input in \"Share Rheumatoid arthritis overview\" dialog")), undefined, "Share it: type a colleague's name, keep View and use");
    await session.step(73, "And user clicks on \"bddsecond\" table row in \"Share Rheumatoid arthritis overview\" dialog", () => clickOn(page, el("\"bddsecond\" table row in \"Share Rheumatoid arthritis overview\" dialog")));
    await session.step(75, "And user unchecks \"Send notifications\" input in \"Share Rheumatoid arthritis overview\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share Rheumatoid arthritis overview\" dialog")), undefined, "Uncheck Send notifications to share without a message");
    await session.step(76, "And user clicks on OK button in \"Share Rheumatoid arthritis overview\" dialog", () => clickOn(page, el("OK button in \"Share Rheumatoid arthritis overview\" dialog")));
    await session.step(77, "Then the \"Share Rheumatoid arthritis overview\" dialog should close", () => dialogCloses(page, "Share Rheumatoid arthritis overview"));
    await session.step(78, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
    await session.step(80, "When user clicks on browse tab", () => clickOn(page, el("browse tab")), undefined, "Done. To find it later: Browse > Dashboards");
    await session.step(81, "Then browse panel should be visible", () => shouldBe(page, el("browse panel"), "visible"));
    await session.step(82, "When user collapses Files tree node inside browse tree", () => collapse(page, el("Files tree node inside browse tree")));
    await session.step(83, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(84, "And user enters \"Rheumatoid arthritis overview\" into gallery search", () => enterInto(page, "Rheumatoid arthritis overview", el("gallery search")));
    await session.step(85, "Then \"Rheumatoid arthritis overview\" gallery card should be visible", () => shouldBe(page, el("\"Rheumatoid arthritis overview\" gallery card"), "visible"));
    await session.step(86, "When the video says \"There it is, for you and for them\"", () => videoSays(page, "There it is, for you and for them"));
    await session.step(87, "And user hovers over \"Rheumatoid arthritis overview\" gallery card", () => hoverOver(page, el("\"Rheumatoid arthritis overview\" gallery card")));
    await session.step(88, "Then no errors should have been logged", () => noErrors(page));
  });
});
