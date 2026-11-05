/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/models/models-gallery.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.models, ml.menu.models.train-model, ml.menu.models.apply-model]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, clickOnHolding, enterInto, fileMatchingDownloaded, isExpanded, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldNotContainText, typeInto, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, columnsExactly} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newestMatchingFilled, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, galleryCountLower, modelsOnServer, noModelOnServer, openTableOf, rememberGalleryCount, switchView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Predictive models gallery", () => {
  const session = feature(test, "features/models/models-gallery.feature", import.meta.url);
  test("The Predictive models gallery", {tag: ["@journey", "@eda", "@realizes:views.models", "@realizes:ml.menu.models.train-model", "@realizes:ml.menu.models.apply-model"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And no predictive model named \"BDD-Kestrel-{run}\" is on the server", () => noModelOnServer(page, session.text("BDD-Kestrel-{run}")));
    await session.step(32, "And no predictive model named \"BDD-Osprey-{run}\" is on the server", () => noModelOnServer(page, session.text("BDD-Osprey-{run}")));
    await run.scenario("Two models are trained and saved, each on a table of its own", async () => {
      await session.step(35, "Given user opens a table \"readings\" with:", () => openTableOf(page, "readings", [["f1","f2","target"],["1","7","A"],["2","3","A"],["3","9","B"],["4","1","A"],["5","8","B"],["6","2","A"],["7","10","B"],["8","4","B"],["9","6","A"],["10","5","B"]]), [["f1","f2","target"],["1","7","A"],["2","3","A"],["3","9","B"],["4","1","A"],["5","8","B"],["6","2","A"],["7","10","B"],["8","4","B"],["9","6","A"],["10","5","B"]]);
      await session.step(47, "When user picks \"ML > Models > Train Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Train Model..."));
      await session.step(48, "And user selects \"target\" in Predict input", () => selectIn(page, "target", el("Predict input")));
      await session.step(49, "And user clicks on editor of Features input", () => clickOn(page, el("editor of Features input")));
      await session.step(50, "And user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(51, "And user clicks on the \"cell 1 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 1 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(52, "And user clicks on the \"cell 2 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 2 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(53, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(54, "Then model preview should be ready", () => shouldBe(page, el("model preview"), "ready"));
      await session.step(55, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(56, "And user enters \"BDD-Kestrel-{run}\" into Name input in dialog", () => enterInto(page, session.text("BDD-Kestrel-{run}"), el("Name input in dialog")));
      await session.step(57, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
      await session.step(58, "Then dialog should be absent", () => shouldBe(page, el("dialog"), "absent"));
      await session.step(59, "And 1 predictive model named \"BDD-Kestrel-{run}\" should be on the server", () => modelsOnServer(page, 1, session.text("BDD-Kestrel-{run}")));
      await session.step(60, "Given user opens a table \"levels\" with:", () => openTableOf(page, "levels", [["g1","g2","score"],["1","1.1","3.1"],["2","2.3","1.2"],["3","2.9","4.4"],["4","4.2","2.5"],["5","5.1","5.3"],["6","6.4","1.7"],["7","6.8","3.9"],["8","8.3","2.2"]]), [["g1","g2","score"],["1","1.1","3.1"],["2","2.3","1.2"],["3","2.9","4.4"],["4","4.2","2.5"],["5","5.1","5.3"],["6","6.4","1.7"],["7","6.8","3.9"],["8","8.3","2.2"]]);
      await session.step(70, "When user picks \"ML > Models > Train Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Train Model..."));
      await session.step(71, "And user selects \"score\" in Predict input", () => selectIn(page, "score", el("Predict input")));
      await session.step(72, "And user clicks on editor of Features input", () => clickOn(page, el("editor of Features input")));
      await session.step(73, "And user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(74, "And user clicks on the \"cell 1 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 1 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(75, "And user clicks on the \"cell 2 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 2 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(76, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(77, "Then model preview should be ready", () => shouldBe(page, el("model preview"), "ready"));
      await session.step(78, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(79, "And user enters \"BDD-Osprey-{run}\" into Name input in dialog", () => enterInto(page, session.text("BDD-Osprey-{run}"), el("Name input in dialog")));
      await session.step(80, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
      await session.step(81, "Then dialog should be absent", () => shouldBe(page, el("dialog"), "absent"));
      await session.step(82, "And 1 predictive model named \"BDD-Osprey-{run}\" should be on the server", () => modelsOnServer(page, 1, session.text("BDD-Osprey-{run}")));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The gallery search finds a model by a word of its name", async () => {
      await session.step(86, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(87, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(88, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
      await session.step(89, "When user clicks on \"Predictive models\" tree node inside browse tree", () => clickOn(page, el("\"Predictive models\" tree node inside browse tree")));
      await session.step(90, "Then the \"Models\" view should be current", () => viewIsCurrent(page, "Models"));
      await session.step(91, "And \"BDD-Kestrel-{run}\" label in gallery should be visible", () => shouldBe(page, el(session.text("\"BDD-Kestrel-{run}\" label in gallery")), "visible"));
      await session.step(92, "And \"BDD-Osprey-{run}\" label in gallery should be visible", () => shouldBe(page, el(session.text("\"BDD-Osprey-{run}\" label in gallery")), "visible"));
      await session.step(93, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(94, "And user types \"Kestrel\" into gallery search", () => typeInto(page, "Kestrel", el("gallery search")));
      await session.step(95, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(96, "And \"BDD-Kestrel-{run}\" label in gallery should be visible", () => shouldBe(page, el(session.text("\"BDD-Kestrel-{run}\" label in gallery")), "visible"));
      await session.step(97, "And \"BDD-Osprey-{run}\" label in gallery should be absent", () => shouldBe(page, el(session.text("\"BDD-Osprey-{run}\" label in gallery")), "absent"));
      await session.step(98, "When user types \"Osprey\" into gallery search", () => typeInto(page, "Osprey", el("gallery search")));
      await session.step(99, "Then \"BDD-Osprey-{run}\" label in gallery should be visible", () => shouldBe(page, el(session.text("\"BDD-Osprey-{run}\" label in gallery")), "visible"));
      await session.step(100, "And the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(101, "And \"BDD-Kestrel-{run}\" label in gallery should be absent", () => shouldBe(page, el(session.text("\"BDD-Kestrel-{run}\" label in gallery")), "absent"));
      await session.step(103, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(104, "Then \"BDD-Kestrel-{run}\" label in gallery should be visible", () => shouldBe(page, el(session.text("\"BDD-Kestrel-{run}\" label in gallery")), "visible"));
      await session.step(105, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Created by me quick filter writes its query into the search", async () => {
      await session.step(108, "When user clicks on \"Toggle filters\" icon", () => clickOn(page, el("\"Toggle filters\" icon")));
      await session.step(109, "Then \"Is applicable to...\" tag should be visible", () => shouldBe(page, el("\"Is applicable to...\" tag"), "visible"));
      await session.step(110, "When user clicks on \"Created by me\" tag", () => clickOn(page, el("\"Created by me\" tag")));
      await session.step(111, "Then gallery search should have value \"author = @current\"", () => shouldHaveValue(page, el("gallery search"), "author = @current"));
      await session.step(112, "When user clicks on \"All\" tag", () => clickOn(page, el("\"All\" tag")));
      await session.step(113, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
      await session.step(114, "When user clicks on \"Toggle filters\" icon", () => clickOn(page, el("\"Toggle filters\" icon")));
      await session.step(115, "Then \"Is applicable to...\" tag should be hidden", () => shouldBe(page, el("\"Is applicable to...\" tag"), "hidden"));
      await session.step(116, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Edit... changes the description the Details pane shows", async () => {
      await session.step(119, "When user picks \"Edit...\" from the context menu of \"BDD-Kestrel-{run}\" label in gallery", () => pickFromContextMenu(page, "Edit...", el(session.text("\"BDD-Kestrel-{run}\" label in gallery"))));
      await session.step(120, "Then \"Predictive model\" dialog should be visible", () => shouldBe(page, el("\"Predictive model\" dialog"), "visible"));
      await session.step(121, "And Name input in \"Predictive model\" dialog should have value \"BDD-Kestrel-{run}\"", () => shouldHaveValue(page, el("Name input in \"Predictive model\" dialog"), session.text("BDD-Kestrel-{run}")));
      await session.step(122, "When user enters \"Classifies the readings, {run}\" into Description input in \"Predictive model\" dialog", () => enterInto(page, session.text("Classifies the readings, {run}"), el("Description input in \"Predictive model\" dialog")));
      await session.step(123, "And user clicks on OK button in \"Predictive model\" dialog", () => clickOn(page, el("OK button in \"Predictive model\" dialog")));
      await session.step(124, "Then the \"Predictive model\" dialog should close", () => dialogCloses(page, "Predictive model"));
      await session.step(125, "When user clicks on \"BDD-Kestrel-{run}\" label in gallery", () => clickOn(page, el(session.text("\"BDD-Kestrel-{run}\" label in gallery"))));
      await session.step(126, "Then the context panel should show \"BDD-Kestrel-{run}\"", () => contextPanelShows(page, session.text("BDD-Kestrel-{run}")));
      await session.step(127, "And \"Details\" pane in context panel should contain text \"Classifies the readings, {run}\"", () => shouldContainText(page, el("\"Details\" pane in context panel"), session.text("Classifies the readings, {run}")));
      await session.step(128, "And \"Details\" pane in context panel should contain text \"f1, f2\"", () => shouldContainText(page, el("\"Details\" pane in context panel"), "f1, f2"));
      await session.step(129, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Run Evaluation draws the model's charts and metrics on its training table", async () => {
      await session.step(132, "When user clicks on \"Performance\" pane in context panel", () => clickOn(page, el("\"Performance\" pane in context panel")));
      await session.step(133, "Then \"Run Evaluation\" button in context panel should be visible", () => shouldBe(page, el("\"Run Evaluation\" button in context panel"), "visible"));
      await session.step(134, "And \"Performance\" pane in context panel should not contain text \"Accuracy\"", () => shouldNotContainText(page, el("\"Performance\" pane in context panel"), "Accuracy"));
      await session.step(135, "When user clicks on \"Run Evaluation\" button in context panel", () => clickOn(page, el("\"Run Evaluation\" button in context panel")));
      await session.step(136, "Then \"Performance\" pane in context panel should contain text \"Accuracy\"", () => shouldContainText(page, el("\"Performance\" pane in context panel"), "Accuracy"));
      await session.step(137, "And \"Performance\" pane in context panel should contain text \"Confusions\"", () => shouldContainText(page, el("\"Performance\" pane in context panel"), "Confusions"));
      await session.step(138, "And scatter plot viewer in \"Performance\" pane in context panel should be visible", () => shouldBe(page, el("scatter plot viewer in \"Performance\" pane in context panel"), "visible"));
      await session.step(139, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(141, "When user clicks on \"Sharing\" pane in context panel", () => clickOn(page, el("\"Sharing\" pane in context panel")));
      await session.step(142, "And user clicks on \"BDD-Osprey-{run}\" gallery card", () => clickOn(page, el(session.text("\"BDD-Osprey-{run}\" gallery card"))));
      await session.step(143, "Then the context panel should show \"BDD-Osprey-{run}\"", () => contextPanelShows(page, session.text("BDD-Osprey-{run}")));
      await session.step(144, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A card applies its model to an open table it fits, and only to that one", async () => {
      await session.step(147, "Then \"BDD-Osprey-{run}\" gallery card should contain text \"Applicable to levels\"", () => shouldContainText(page, el(session.text("\"BDD-Osprey-{run}\" gallery card")), "Applicable to levels"));
      await session.step(148, "And \"BDD-Osprey-{run}\" gallery card should not contain text \"readings\"", () => shouldNotContainText(page, el(session.text("\"BDD-Osprey-{run}\" gallery card")), "readings"));
      await session.step(149, "And \"BDD-Kestrel-{run}\" gallery card should contain text \"Applicable to readings\"", () => shouldContainText(page, el(session.text("\"BDD-Kestrel-{run}\" gallery card")), "Applicable to readings"));
      await session.step(150, "And \"BDD-Kestrel-{run}\" gallery card should not contain text \"levels\"", () => shouldNotContainText(page, el(session.text("\"BDD-Kestrel-{run}\" gallery card")), "levels"));
      await session.step(151, "When user picks \"Apply to > levels (8 rows, 3 columns)\" from the context menu of \"BDD-Osprey-{run}\" label in gallery", () => pickFromContextMenu(page, "Apply to > levels (8 rows, 3 columns)", el(session.text("\"BDD-Osprey-{run}\" label in gallery"))));
      await session.step(152, "Then \"Apply predictive model\" dialog should be visible", () => shouldBe(page, el("\"Apply predictive model\" dialog"), "visible"));
      await session.step(153, "And Model input in \"Apply predictive model\" dialog should contain text \"BDD-Osprey-\"", () => shouldContainText(page, el("Model input in \"Apply predictive model\" dialog"), "BDD-Osprey-"));
      await session.step(154, "And Model input in \"Apply predictive model\" dialog should not contain text \"BDD-Kestrel-\"", () => shouldNotContainText(page, el("Model input in \"Apply predictive model\" dialog"), "BDD-Kestrel-"));
      await session.step(155, "And Inputs input in \"Apply predictive model\" dialog should contain text \"(2/2)\"", () => shouldContainText(page, el("Inputs input in \"Apply predictive model\" dialog"), "(2/2)"));
      await session.step(156, "When user clicks on OK button in \"Apply predictive model\" dialog", () => clickOn(page, el("OK button in \"Apply predictive model\" dialog")));
      await session.step(157, "Then the \"Apply predictive model\" dialog should close", () => dialogCloses(page, "Apply predictive model"));
      await session.step(158, "And the \"levels\" view should be current", () => viewIsCurrent(page, "levels"));
      await session.step(159, "And the table should have 4 columns", () => columnCount(page, 4));
      await session.step(160, "And the newest column matching \"^score \\(\" should have no missing values", () => newestMatchingFilled(page, "^score \\("));
      await session.step(161, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(162, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Save as Zip downloads the model", async () => {
      await session.step(165, "Given user watches downloads", () => watchDownloads(page));
      await session.step(166, "And user switches to the \"Models\" view", () => switchView(page, "Models"));
      await session.step(167, "When user picks \"Save as Zip\" from the context menu of \"BDD-Kestrel-{run}\" label in gallery", () => pickFromContextMenu(page, "Save as Zip", el(session.text("\"BDD-Kestrel-{run}\" label in gallery"))));
      await session.step(168, "Then a file matching \"\\.zip$\" should have been downloaded", () => fileMatchingDownloaded(page, "\\.zip$"));
      await session.step(169, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Two cards picked with Ctrl+click are compared side by side", async () => {
      await session.step(172, "When user clicks on \"BDD-Osprey-{run}\" gallery card", () => clickOn(page, el(session.text("\"BDD-Osprey-{run}\" gallery card"))));
      await session.step(173, "Then the context panel should show \"BDD-Osprey-{run}\"", () => contextPanelShows(page, session.text("BDD-Osprey-{run}")));
      await session.step(174, "When user clicks on \"BDD-Kestrel-{run}\" gallery card holding Control", () => clickOnHolding(page, el(session.text("\"BDD-Kestrel-{run}\" gallery card")), "Control"));
      await session.step(175, "Then context panel should contain text \"2 models\"", () => shouldContainText(page, el("context panel"), "2 models"));
      await session.step(176, "When user clicks on \"Actions\" pane in context panel", () => clickOn(page, el("\"Actions\" pane in context panel")));
      await session.step(177, "And user clicks on \"Compare\" label in \"Actions\" pane in context panel", () => clickOn(page, el("\"Compare\" label in \"Actions\" pane in context panel")));
      await session.step(178, "Then the \"Compare models\" view should be current", () => viewIsCurrent(page, "Compare models"));
      await session.step(179, "And the table should have 2 rows", () => rowCount(page, 2));
      await session.step(180, "And the table should have the columns \"Name, Description, Method, Source\"", () => columnsExactly(page, "Name, Description, Method, Source"));
      await session.step(181, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
