/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/eda/dashboards.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.dashboards]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {ownConnectionGone, ownProjectGone, startTutorial, stepDone, stepDoneTimes, stepNotDone, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, expand, insertLine} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, noHintShown, projectsOnServer, queriesOnServer, standReachesConnection, standRunsService, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Dashboards tutorial", () => {
  const session = feature(test, "features/eda/dashboards.feature", import.meta.url);
  test("A learner completes the Dashboards tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.dashboards"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(27, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(28, "And the user's own project \"Coffee sales dashboard\" is removed now and at feature end", () => ownProjectGone(page, "Coffee sales dashboard"));
    await session.step(29, "And the user's own connection \"Starbucks\" and query \"Stores in @state\" are removed now and at feature end", () => ownConnectionGone(page, "Starbucks", "Stores in @state"));
    await session.step(30, "And the \"Dashboards\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Dashboards"));
    await session.step(31, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(35, "Given the stand runs the \"Grok Connect\" service", () => standRunsService(page, "Grok Connect"));
    await session.step(36, "When user starts the \"Dashboards\" tutorial", () => startTutorial(page, "Dashboards"));
    await session.step(37, "Then the tutorial progress should be 1 of 27", () => tutorialProgress(page, 1, 27));
    await session.step(38, "Given the tutorial step \"Create a connection to Postgres server\" should not be done yet", () => stepNotDone(page, "Create a connection to Postgres server"));
    await session.step(39, "When user picks \"New connection...\" from the context menu of Databases---Postgres tree node inside browse tree", () => pickFromContextMenu(page, "New connection...", el("Databases---Postgres tree node inside browse tree")));
    await session.step(40, "Then the tutorial step \"Create a connection to Postgres server\" should be done", () => stepDone(page, "Create a connection to Postgres server"));
    await session.step(41, "When user enters \"Starbucks\" into \"Name\" input in \"Add new connection\" dialog", () => enterInto(page, "Starbucks", el("\"Name\" input in \"Add new connection\" dialog")));
    await session.step(42, "Then the tutorial step \"Set \\\"Name\\\" to \\\"Starbucks\\\"\" should be done", () => stepDone(page, "Set \"Name\" to \"Starbucks\""));
    await session.step(43, "When user enters \"db.datagrok.ai\" into \"Server\" input in \"Add new connection\" dialog", () => enterInto(page, "db.datagrok.ai", el("\"Server\" input in \"Add new connection\" dialog")));
    await session.step(44, "And user enters \"54324\" into \"Port\" input in \"Add new connection\" dialog", () => enterInto(page, "54324", el("\"Port\" input in \"Add new connection\" dialog")));
    await session.step(45, "And user enters \"starbucks\" into \"Db\" input in \"Add new connection\" dialog", () => enterInto(page, "starbucks", el("\"Db\" input in \"Add new connection\" dialog")));
    await session.step(46, "And user enters \"datagrok\" into \"Login\" input in \"Add new connection\" dialog", () => enterInto(page, "datagrok", el("\"Login\" input in \"Add new connection\" dialog")));
    await session.step(47, "And user enters \"KKfIh6ooS7vjzHYrNiRrderyz3KUyglrhSJF\" into \"Password\" input in \"Add new connection\" dialog", () => enterInto(page, "KKfIh6ooS7vjzHYrNiRrderyz3KUyglrhSJF", el("\"Password\" input in \"Add new connection\" dialog")));
    await session.step(48, "Then the tutorial step \"Set \\\"Password\\\" to \\\"KKfIh6ooS7vjzHYrNiRrderyz3KUyglrhSJF\\\"\" should be done", () => stepDone(page, "Set \"Password\" to \"KKfIh6ooS7vjzHYrNiRrderyz3KUyglrhSJF\""));
    await session.step(49, "When user clicks on OK button in \"Add new connection\" dialog", () => clickOn(page, el("OK button in \"Add new connection\" dialog")));
    await session.step(50, "Then the tutorial step \"Click \\\"OK\\\"\" should be done", () => stepDone(page, "Click \"OK\""));
    await session.step(51, "Given the tutorial step \"Create a data query to the \\\"Starbucks\\\" data connection\" should not be done yet", () => stepNotDone(page, "Create a data query to the \"Starbucks\" data connection"));
    await session.step(52, "When user picks \"New Query...\" from the context menu of Databases---Postgres---Starbucks tree node inside browse tree", () => pickFromContextMenu(page, "New Query...", el("Databases---Postgres---Starbucks tree node inside browse tree")));
    await session.step(53, "Then the tutorial step \"Create a data query to the \\\"Starbucks\\\" data connection\" should be done", () => stepDone(page, "Create a data query to the \"Starbucks\" data connection"));
    await session.step(54, "Given the tutorial step \"Set \\\"Name\\\" to \\\"Stores in @state\\\"\" should not be done yet", () => stepNotDone(page, "Set \"Name\" to \"Stores in @state\""));
    await session.step(55, "When user enters \"Stores in @state\" into \"Name\" input", () => enterInto(page, "Stores in @state", el("\"Name\" input")));
    await session.step(56, "Then the tutorial step \"Set \\\"Name\\\" to \\\"Stores in @state\\\"\" should be done", () => stepDone(page, "Set \"Name\" to \"Stores in @state\""));
    await session.step(57, "When user puts \"select * from starbucks_us where state = @state;\" on the first line of code editor", () => insertLine(page, "select * from starbucks_us where state = @state;", el("code editor")));
    await session.step(58, "Then the tutorial step \"Add \\\"select * from starbucks_us where state = @state;\\\" to the editor\" should be done", () => stepDone(page, "Add \"select * from starbucks_us where state = @state;\" to the editor"));
    await session.step(59, "When user puts \"--input: string state\" on the first line of code editor", () => insertLine(page, "--input: string state", el("code editor")));
    await session.step(60, "Then the tutorial step \"Add \\\"--input: string state\\\" as the first line of the query\" should be done", () => stepDone(page, "Add \"--input: string state\" as the first line of the query"));
    await session.step(61, "When user clicks on SAVE button", () => clickOn(page, el("SAVE button")));
    await session.step(62, "Then the tutorial step \"Save the query\" should be done", () => stepDone(page, "Save the query"));
    await session.step(63, "And 1 query named \"Stores in @state\" should be on the server", () => queriesOnServer(page, 1, "Stores in @state"));
    await session.step(64, "Given the tutorial step \"Find Browse on the sidebar and click\" should not be done yet", () => stepNotDone(page, "Find Browse on the sidebar and click"));
    await session.step(65, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(66, "Then the tutorial step \"Find Browse on the sidebar and click\" should be done", () => stepDone(page, "Find Browse on the sidebar and click"));
    await session.step(67, "Given the tutorial step \"Find the created query in the browse view, right-click it and hit Run\" should not be done yet", () => stepNotDone(page, "Find the created query in the browse view, right-click it and hit Run"));
    await session.step(68, "When user expands Databases---Postgres---Starbucks tree node inside browse tree", () => expand(page, el("Databases---Postgres---Starbucks tree node inside browse tree")));
    await session.step(69, "And user picks \"Run\" from the context menu of Databases---Postgres---Starbucks---Stores-in-@state tree node inside browse tree", () => pickFromContextMenu(page, "Run", el("Databases---Postgres---Starbucks---Stores-in-@state tree node inside browse tree")));
    await session.step(70, "Then the tutorial step \"Find the created query in the browse view, right-click it and hit Run\" should be done", () => stepDone(page, "Find the created query in the browse view, right-click it and hit Run"));
    await session.step(71, "When user enters \"NY\" into \"State\" input in \"Stores in @state\" dialog", () => enterInto(page, "NY", el("\"State\" input in \"Stores in @state\" dialog")));
    await session.step(72, "Then the tutorial step \"Set state to \\\"NY\\\"\" should be done", () => stepDone(page, "Set state to \"NY\""));
    await session.step(75, "Given the stand can reach the database of the \"Starbucks\" connection", () => standReachesConnection(page, "Starbucks"));
    await session.step(76, "When user clicks on OK button in \"Stores in @state\" dialog", () => clickOn(page, el("OK button in \"Stores in @state\" dialog")));
    await session.step(77, "Then the tutorial step \"Click \\\"OK\\\" to run the query\" should be done", () => stepDone(page, "Click \"OK\" to run the query"));
    await session.step(78, "And the table should have 645 rows", () => rowCount(page, 645));
    await session.step(80, "Given the tutorial step \"Open bar chart\" should not be done yet", () => stepNotDone(page, "Open bar chart"));
    await session.step(81, "When user clicks on bar-chart icon in toolbox", () => clickOn(page, el("bar-chart icon in toolbox")));
    await session.step(82, "Then the tutorial step \"Open bar chart\" should be done", () => stepDone(page, "Open bar chart"));
    await session.step(83, "Given the tutorial step \"Save a project\" should not be done yet", () => stepNotDone(page, "Save a project"));
    await session.step(84, "When user clicks on SAVE button", () => clickOn(page, el("SAVE button")));
    await session.step(85, "Then the tutorial step \"Save a project\" should be done", () => stepDone(page, "Save a project"));
    await session.step(86, "When user enters \"Coffee sales dashboard\" into \"Name\" text input in \"Save project\" dialog", () => enterInto(page, "Coffee sales dashboard", el("\"Name\" text input in \"Save project\" dialog")));
    await session.step(87, "Then the tutorial step \"Set the project name to \\\"Coffee sales dashboard\\\"\" should be done", () => stepDone(page, "Set the project name to \"Coffee sales dashboard\""));
    await session.step(88, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(89, "Then the tutorial step \"Click \\\"OK\\\"\" should be done 2 times", () => stepDoneTimes(page, "Click \"OK\"", 2));
    await session.step(90, "When user clicks on CANCEL button in \"Share Coffee sales dashboard\" dialog", () => clickOn(page, el("CANCEL button in \"Share Coffee sales dashboard\" dialog")));
    await session.step(91, "Then the tutorial step \"Skip the sharing step\" should be done", () => stepDone(page, "Skip the sharing step"));
    await session.step(92, "And 1 project named \"Coffee sales dashboard\" should be on the server", () => projectsOnServer(page, 1, "Coffee sales dashboard"));
    await session.step(93, "Given the tutorial step \"Close the project\" should not be done yet", () => stepNotDone(page, "Close the project"));
    await session.step(94, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
    await session.step(95, "Then the tutorial step \"Close the project\" should be done", () => stepDone(page, "Close the project"));
    await session.step(96, "Given the tutorial step \"Open browse and click on Dashboards\" should not be done yet", () => stepNotDone(page, "Open browse and click on Dashboards"));
    await session.step(97, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(98, "Then the tutorial step \"Open browse and click on Dashboards\" should be done", () => stepDone(page, "Open browse and click on Dashboards"));
    await session.step(99, "When user double-clicks on \"Coffee sales dashboard\" gallery card", () => doubleClickOn(page, el("\"Coffee sales dashboard\" gallery card")));
    await session.step(100, "Then the tutorial step \"Find and open your project\" should be done", () => stepDone(page, "Find and open your project"));
    await session.step(101, "Given the tutorial step \"Set State to LA\" should not be done yet", () => stepNotDone(page, "Set State to LA"));
    await session.step(102, "When user enters \"LA\" into \"State\" input", () => enterInto(page, "LA", el("\"State\" input")));
    await session.step(103, "Then the tutorial step \"Set State to LA\" should be done", () => stepDone(page, "Set State to LA"));
    await session.step(104, "When user clicks on REFRESH button", () => clickOn(page, el("REFRESH button")));
    await session.step(105, "Then the tutorial step \"Click REFRESH button\" should be done", () => stepDone(page, "Click REFRESH button"));
    await session.step(107, "And the table should have 84 rows", () => rowCount(page, 84));
    await session.step(109, "And the \"Dashboards\" tutorial should be completed", () => tutorialCompleted(page, "Dashboards"));
    await session.step(110, "And the tutorial should have listed 27 steps", () => tutorialStepsListed(page, 27));
    await session.step(111, "And the tutorial progress should be 27 of 27", () => tutorialProgress(page, 27, 27));
    await session.step(112, "And no hint should be shown", () => noHintShown(page));
    await session.step(113, "And no errors should have been logged", () => noErrors(page));
  });
});
