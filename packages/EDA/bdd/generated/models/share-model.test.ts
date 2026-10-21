/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/models/share-model.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sharing.share-dialog, sharing.browse-shared-with-me]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, enterInto, isExpanded, selectIn, shouldBe, shouldContainText, shouldNotContainText, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {newColumnMatching, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, expandShared, modelsOnServer, noModelOnServer, openDataset, pickSharingUser, refreshBrowse, removeSharingUser, runningAccountSignedIn, selectModel, sharedShouldBe, sharingPaneLists, sharingPaneListsNot, sharingUserLevel, sharingUserSignedIn, signInAsSelf, signInAsSharingUser, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, pickFromContextMenu, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Sharing a predictive model", () => {
  const session = feature(test, "features/models/share-model.feature", import.meta.url);
  test("Sharing a predictive model", {tag: ["@journey", "@serial", "@eda", "@realizes:sharing.share-dialog", "@realizes:sharing.browse-shared-with-me"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And no predictive model named \"BDD-Share-Model-{run}\" is on the server", () => noModelOnServer(page, session.text("BDD-Share-Model-{run}")));
    await run.scenario("The owner trains a model and saves it", async () => {
      await session.step(34, "Given user opens demog dataset", () => openDataset(page, ds("demog")));
      await session.step(35, "When user picks \"ML > Models > Train Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Train Model..."));
      await session.step(36, "And user selects \"SEX\" in Predict input", () => selectIn(page, "SEX", el("Predict input")));
      await session.step(37, "And user clicks on editor of Features input", () => clickOn(page, el("editor of Features input")));
      await session.step(38, "And user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(39, "Then the \"text of cell 6 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"HEIGHT\"", () => readingReads(page, "text of cell 6 of __name", el("grid viewer in \"Select columns...\" dialog"), "HEIGHT"));
      await session.step(40, "And the \"text of cell 7 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"WEIGHT\"", () => readingReads(page, "text of cell 7 of __name", el("grid viewer in \"Select columns...\" dialog"), "WEIGHT"));
      await session.step(41, "When user clicks on the \"cell 6 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 6 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(42, "And user clicks on the \"cell 7 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 7 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(43, "Then \"Select columns...\" dialog should contain text \"2 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "2 checked"));
      await session.step(44, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(45, "And user checks \"Ignore missing\" input", () => check(page, el("\"Ignore missing\" input")));
      await session.step(46, "And user selects \"Eda: XGBoost\" in \"Model Engine\" input", () => selectIn(page, "Eda: XGBoost", el("\"Model Engine\" input")));
      await session.step(47, "Then model preview should be ready", () => shouldBe(page, el("model preview"), "ready"));
      await session.step(48, "And \"Accuracy\" table row should be visible", () => shouldBe(page, el("\"Accuracy\" table row"), "visible"));
      await session.step(49, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(50, "And user enters \"BDD-Share-Model-{run}\" into Name input in dialog", () => enterInto(page, session.text("BDD-Share-Model-{run}"), el("Name input in dialog")));
      await session.step(51, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
      await session.step(52, "Then 1 predictive model named \"BDD-Share-Model-{run}\" should be on the server", () => modelsOnServer(page, 1, session.text("BDD-Share-Model-{run}")));
    });
    await run.scenario("The saved model is the owner's alone", async () => {
      await session.step(55, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(56, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(57, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
      await session.step(58, "When user clicks on \"Predictive models\" tree node inside browse tree", () => clickOn(page, el("\"Predictive models\" tree node inside browse tree")));
      await session.step(59, "Then the \"Models\" view should be current", () => viewIsCurrent(page, "Models"));
      await session.step(60, "When user clicks on \"BDD-Share-Model-{run}\" label in gallery", () => clickOn(page, el(session.text("\"BDD-Share-Model-{run}\" label in gallery"))));
      await session.step(61, "Then the context panel should show \"BDD-Share-Model-{run}\"", () => contextPanelShows(page, session.text("BDD-Share-Model-{run}")));
      await session.step(62, "And the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
    });
    await run.scenario("The Share dialog asks who, how much, and whether to notify", async () => {
      await session.step(65, "When user picks \"Share...\" from the context menu of \"BDD-Share-Model-{run}\" label in gallery", () => pickFromContextMenu(page, "Share...", el(session.text("\"BDD-Share-Model-{run}\" label in gallery"))));
      await session.step(66, "Then \"Share BDD-Share-Model-{run}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDD-Share-Model-{run}\" dialog")), "visible"));
      await session.step(67, "And \"User, group, or email\" input in \"Share BDD-Share-Model-{run}\" dialog should be visible", () => shouldBe(page, el(session.text("\"User, group, or email\" input in \"Share BDD-Share-Model-{run}\" dialog")), "visible"));
      await session.step(68, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(69, "And \"Send notifications\" input in \"Share BDD-Share-Model-{run}\" dialog should be hidden", () => shouldBe(page, el(session.text("\"Send notifications\" input in \"Share BDD-Share-Model-{run}\" dialog")), "hidden"));
      await session.step(70, "And \"Share BDD-Share-Model-{run}\" dialog should not contain text \"will also be shared\"", () => shouldNotContainText(page, el(session.text("\"Share BDD-Share-Model-{run}\" dialog")), "will also be shared"));
      await session.step(71, "When user clicks on CANCEL button in \"Share BDD-Share-Model-{run}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDD-Share-Model-{run}\" dialog"))));
      await session.step(72, "Then \"Share BDD-Share-Model-{run}\" dialog should be hidden", () => shouldBe(page, el(session.text("\"Share BDD-Share-Model-{run}\" dialog")), "hidden"));
      await session.step(73, "When user clicks on \"BDD-Share-Model-{run}\" label in gallery", () => clickOn(page, el(session.text("\"BDD-Share-Model-{run}\" label in gallery"))));
      await session.step(74, "Then the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
    });
    await run.scenario("The model is shared with the second account to view and use", async () => {
      await session.step(77, "When user picks \"Share...\" from the context menu of \"BDD-Share-Model-{run}\" label in gallery", () => pickFromContextMenu(page, "Share...", el(session.text("\"BDD-Share-Model-{run}\" label in gallery"))));
      await session.step(78, "Then \"Share BDD-Share-Model-{run}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDD-Share-Model-{run}\" dialog")), "visible"));
      await session.step(79, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(80, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDD-Share-Model-{run}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDD-Share-Model-{run}\" dialog"))));
      await session.step(81, "Then \"Send notifications\" input in \"Share BDD-Share-Model-{run}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Send notifications\" input in \"Share BDD-Share-Model-{run}\" dialog")), "visible"));
      await session.step(82, "And \"Share BDD-Share-Model-{run}\" dialog should not contain text \"will also be shared\"", () => shouldNotContainText(page, el(session.text("\"Share BDD-Share-Model-{run}\" dialog")), "will also be shared"));
      await session.step(83, "When user unchecks \"Send notifications\" input in \"Share BDD-Share-Model-{run}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDD-Share-Model-{run}\" dialog"))));
      await session.step(84, "And user clicks on OK button in \"Share BDD-Share-Model-{run}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDD-Share-Model-{run}\" dialog"))));
      await session.step(85, "Then the \"Share BDD-Share-Model-{run}\" dialog should close", () => dialogCloses(page, session.text("Share BDD-Share-Model-{run}")));
      await session.step(86, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(87, "When user clicks on \"BDD-Share-Model-{run}\" label in gallery", () => clickOn(page, el(session.text("\"BDD-Share-Model-{run}\" label in gallery"))));
      await session.step(88, "Then the sharing pane should list the sharing user", () => sharingPaneLists(page));
    });
    await run.scenario("The second account finds the model under Shared with me and applies it", async () => {
      await session.step(91, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(92, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(93, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(94, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(95, "Given \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(96, "And \"My stuff > Shared with me\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > Shared with me\" tree node inside browse tree")));
      await session.step(97, "When user expands \".\" shared by the running account", () => expandShared(page, "."));
      await session.step(98, "Then \"BDD-Share-Model-{run}\" shared by the running account should be visible", () => sharedShouldBe(page, session.text("BDD-Share-Model-{run}"), "visible"));
      await session.step(99, "Given user opens demog dataset", () => openDataset(page, ds("demog")));
      await session.step(100, "When user picks \"ML > Models > Apply Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Apply Model..."));
      await session.step(101, "Then \"Apply predictive model\" dialog should be visible", () => shouldBe(page, el("\"Apply predictive model\" dialog"), "visible"));
      await session.step(102, "When user selects the predictive model \"BDD-Share-Model-{run}\" in Model input in \"Apply predictive model\" dialog", () => selectModel(page, session.text("BDD-Share-Model-{run}"), el("Model input in \"Apply predictive model\" dialog")));
      await session.step(103, "And user clicks on OK button in \"Apply predictive model\" dialog", () => clickOn(page, el("OK button in \"Apply predictive model\" dialog")));
      await session.step(104, "Then the \"Apply predictive model\" dialog should close", () => dialogCloses(page, "Apply predictive model"));
      await session.step(105, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(106, "And a new column matching \"^SEX\" should have been added", () => newColumnMatching(page, "^SEX"));
      await session.step(107, "And no errors should have been logged", () => noErrors(page));
      await session.step(108, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(109, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(110, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(111, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(112, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
      await session.step(113, "When user clicks on \"Predictive models\" tree node inside browse tree", () => clickOn(page, el("\"Predictive models\" tree node inside browse tree")));
      await session.step(114, "Then the \"Models\" view should be current", () => viewIsCurrent(page, "Models"));
    });
    await run.scenario("The owner takes the share back", async () => {
      await session.step(117, "When user picks \"Share...\" from the context menu of \"BDD-Share-Model-{run}\" label in gallery", () => pickFromContextMenu(page, "Share...", el(session.text("\"BDD-Share-Model-{run}\" label in gallery"))));
      await session.step(118, "Then the access level of the sharing user in \"Share BDD-Share-Model-{run}\" dialog should be \"View and use\"", () => sharingUserLevel(page, el(session.text("\"Share BDD-Share-Model-{run}\" dialog")), "View and use"));
      await session.step(119, "When user removes the sharing user from \"Share BDD-Share-Model-{run}\" dialog", () => removeSharingUser(page, el(session.text("\"Share BDD-Share-Model-{run}\" dialog"))));
      await session.step(120, "And user clicks on OK button in \"Share BDD-Share-Model-{run}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDD-Share-Model-{run}\" dialog"))));
      await session.step(121, "Then the \"Share BDD-Share-Model-{run}\" dialog should close", () => dialogCloses(page, session.text("Share BDD-Share-Model-{run}")));
      await session.step(122, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(123, "When user clicks on \"BDD-Share-Model-{run}\" label in gallery", () => clickOn(page, el(session.text("\"BDD-Share-Model-{run}\" label in gallery"))));
      await session.step(124, "Then the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(125, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
