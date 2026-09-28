/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/data-access/data-connectors.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.data-connectors]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {ownConnectionGone, startTutorial, stepDone, stepNotDone, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, insertLine, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, connectionsOnServer, noHintShown, standReachesConnection, standRunsService, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Data Connectors tutorial", () => {
  const session = feature(test, "features/data-access/data-connectors.feature", import.meta.url);
  test("A learner completes the Data Connectors tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.data-connectors"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(28, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(29, "And the user's own connection \"Starbucks\" and query \"Get Starbucks US\" are removed now and at feature end", () => ownConnectionGone(page, "Starbucks", "Get Starbucks US"));
    await session.step(30, "And the \"Data Connectors\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Data Connectors"));
    await session.step(31, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(35, "Given the stand runs the \"Grok Connect\" service", () => standRunsService(page, "Grok Connect"));
    await session.step(36, "When user starts the \"Data Connectors\" tutorial", () => startTutorial(page, "Data Connectors"));
    await session.step(37, "Then the tutorial progress should be 1 of 11", () => tutorialProgress(page, 1, 11));
    await session.step(38, "Given the tutorial step \"Create a connection to Postgres server\" should not be done yet", () => stepNotDone(page, "Create a connection to Postgres server"));
    await session.step(39, "When user picks \"New connection...\" from the context menu of Databases---Postgres tree node inside browse tree", () => pickFromContextMenu(page, "New connection...", el("Databases---Postgres tree node inside browse tree")));
    await session.step(40, "Then the tutorial step \"Create a connection to Postgres server\" should be done", () => stepDone(page, "Create a connection to Postgres server"));
    await session.step(41, "And \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(42, "When user enters \"Starbucks\" into \"Name\" input in \"Add new connection\" dialog", () => enterInto(page, "Starbucks", el("\"Name\" input in \"Add new connection\" dialog")));
    await session.step(43, "Then the tutorial step \"Set \\\"Name\\\" to \\\"Starbucks\\\"\" should be done", () => stepDone(page, "Set \"Name\" to \"Starbucks\""));
    await session.step(44, "When user enters \"db.datagrok.ai\" into \"Server\" input in \"Add new connection\" dialog", () => enterInto(page, "db.datagrok.ai", el("\"Server\" input in \"Add new connection\" dialog")));
    await session.step(45, "Then the tutorial step \"Set \\\"Server\\\" to \\\"db.datagrok.ai\\\"\" should be done", () => stepDone(page, "Set \"Server\" to \"db.datagrok.ai\""));
    await session.step(46, "When user enters \"54324\" into \"Port\" input in \"Add new connection\" dialog", () => enterInto(page, "54324", el("\"Port\" input in \"Add new connection\" dialog")));
    await session.step(47, "Then the tutorial step \"Set \\\"Port\\\" to \\\"54324\\\"\" should be done", () => stepDone(page, "Set \"Port\" to \"54324\""));
    await session.step(48, "When user enters \"starbucks\" into \"Db\" input in \"Add new connection\" dialog", () => enterInto(page, "starbucks", el("\"Db\" input in \"Add new connection\" dialog")));
    await session.step(49, "Then the tutorial step \"Set \\\"Db\\\" to \\\"starbucks\\\"\" should be done", () => stepDone(page, "Set \"Db\" to \"starbucks\""));
    await session.step(50, "When user enters \"datagrok\" into \"Login\" input in \"Add new connection\" dialog", () => enterInto(page, "datagrok", el("\"Login\" input in \"Add new connection\" dialog")));
    await session.step(51, "Then the tutorial step \"Set \\\"Login\\\" to \\\"datagrok\\\"\" should be done", () => stepDone(page, "Set \"Login\" to \"datagrok\""));
    await session.step(52, "When user enters \"KKfIh6ooS7vjzHYrNiRrderyz3KUyglrhSJF\" into \"Password\" input in \"Add new connection\" dialog", () => enterInto(page, "KKfIh6ooS7vjzHYrNiRrderyz3KUyglrhSJF", el("\"Password\" input in \"Add new connection\" dialog")));
    await session.step(53, "Then the tutorial step \"Set \\\"Password\\\" to \\\"KKfIh6ooS7vjzHYrNiRrderyz3KUyglrhSJF\\\"\" should be done", () => stepDone(page, "Set \"Password\" to \"KKfIh6ooS7vjzHYrNiRrderyz3KUyglrhSJF\""));
    await session.step(54, "When user clicks on OK button in \"Add new connection\" dialog", () => clickOn(page, el("OK button in \"Add new connection\" dialog")));
    await session.step(55, "Then the tutorial step \"Click \\\"OK\\\"\" should be done", () => stepDone(page, "Click \"OK\""));
    await session.step(56, "And 1 connection named \"Starbucks\" should be on the server", () => connectionsOnServer(page, 1, "Starbucks"));
    await session.step(58, "Given the tutorial step \"Create a data query to the \\\"Starbucks\\\" data connection\" should not be done yet", () => stepNotDone(page, "Create a data query to the \"Starbucks\" data connection"));
    await session.step(59, "When user picks \"New Query...\" from the context menu of Databases---Postgres---Starbucks tree node inside browse tree", () => pickFromContextMenu(page, "New Query...", el("Databases---Postgres---Starbucks tree node inside browse tree")));
    await session.step(60, "Then the tutorial step \"Create a data query to the \\\"Starbucks\\\" data connection\" should be done", () => stepDone(page, "Create a data query to the \"Starbucks\" data connection"));
    await session.step(61, "Given the tutorial step \"Set \\\"Name\\\" to \\\"Get Starbucks US\\\"\" should not be done yet", () => stepNotDone(page, "Set \"Name\" to \"Get Starbucks US\""));
    await session.step(62, "When user enters \"Get Starbucks US\" into \"Name\" input", () => enterInto(page, "Get Starbucks US", el("\"Name\" input")));
    await session.step(63, "Then the tutorial step \"Set \\\"Name\\\" to \\\"Get Starbucks US\\\"\" should be done", () => stepDone(page, "Set \"Name\" to \"Get Starbucks US\""));
    await session.step(66, "Given the stand can reach the database of the \"Starbucks\" connection", () => standReachesConnection(page, "Starbucks"));
    await session.step(67, "When user puts \"select * from starbucks_us\" on the first line of code editor", () => insertLine(page, "select * from starbucks_us", el("code editor")));
    await session.step(68, "And user clicks on play icon", () => clickOn(page, el("play icon")));
    await session.step(69, "Then the tutorial step \"Add \\\"select * from starbucks_us\\\" to the editor and hit \\\"Play\\\"\" should be done", () => stepDone(page, "Add \"select * from starbucks_us\" to the editor and hit \"Play\""));
    await session.step(70, "And the \"Data Connectors\" tutorial should be completed", () => tutorialCompleted(page, "Data Connectors"));
    await session.step(71, "And the tutorial should have listed 11 steps", () => tutorialStepsListed(page, 11));
    await session.step(72, "And the tutorial progress should be 11 of 11", () => tutorialProgress(page, 11, 11));
    await session.step(73, "And no hint should be shown", () => noHintShown(page));
    await session.step(74, "And no errors should have been logged", () => noErrors(page));
  });
});
