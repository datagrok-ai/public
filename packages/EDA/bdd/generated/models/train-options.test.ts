/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/models/train-options.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [ml.menu.models.train-model, ml.menu.models.apply-model, eda.model.linear-regression]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, enterInto, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, displayedInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {dialogCloses, modelsOnServer, noModelOnServer, openDataset, openTableOf, selectModel, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The preprocessing options of Train Model, kept by the saved model", () => {
  const session = feature(test, "features/models/train-options.feature", import.meta.url);
  test("The preprocessing options of Train Model, kept by the saved model", {tag: ["@journey", "@eda", "@realizes:ml.menu.models.train-model", "@realizes:ml.menu.models.apply-model", "@realizes:eda.model.linear-regression"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And no predictive model named \"BDD-Probability-{run}\" is on the server", () => noModelOnServer(page, session.text("BDD-Probability-{run}")));
    await session.step(27, "And no predictive model named \"BDD-OneHot-{run}\" is on the server", () => noModelOnServer(page, session.text("BDD-OneHot-{run}")));
    await run.scenario("Ignore missing is offered for missing values, and hides Impute missing once ticked", async () => {
      await session.step(30, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(31, "When user picks \"ML > Models > Train Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Train Model..."));
      await session.step(32, "Then the \"Predictive model\" view should be current", () => viewIsCurrent(page, "Predictive model"));
      await session.step(33, "When user selects \"SEX\" in Predict input", () => selectIn(page, "SEX", el("Predict input")));
      await session.step(34, "And user clicks on editor of Features input", () => clickOn(page, el("editor of Features input")));
      await session.step(35, "And user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(36, "Then the \"text of cell 6 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"HEIGHT\"", () => readingReads(page, "text of cell 6 of __name", el("grid viewer in \"Select columns...\" dialog"), "HEIGHT"));
      await session.step(37, "And the \"text of cell 7 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"WEIGHT\"", () => readingReads(page, "text of cell 7 of __name", el("grid viewer in \"Select columns...\" dialog"), "WEIGHT"));
      await session.step(38, "When user clicks on the \"cell 6 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 6 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(39, "And user clicks on the \"cell 7 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 7 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(40, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(41, "Then editor of Features input should contain text \"(2)\"", () => shouldContainText(page, el("editor of Features input"), "(2)"));
      await session.step(42, "And model preview should contain text \"Column 'HEIGHT' contains missing values.\"", () => shouldContainText(page, el("model preview"), "Column 'HEIGHT' contains missing values."));
      await session.step(43, "And model preview should be invalid", () => shouldBe(page, el("model preview"), "invalid"));
      await session.step(44, "And \"Ignore missing\" input should be visible", () => shouldBe(page, el("\"Ignore missing\" input"), "visible"));
      await session.step(45, "And \"Impute missing\" input should be visible", () => shouldBe(page, el("\"Impute missing\" input"), "visible"));
      await session.step(46, "And \"Predict probability\" input should be visible", () => shouldBe(page, el("\"Predict probability\" input"), "visible"));
      await session.step(47, "And \"Model Engine\" input should be absent", () => shouldBe(page, el("\"Model Engine\" input"), "absent"));
      await session.step(48, "When user checks \"Ignore missing\" input", () => check(page, el("\"Ignore missing\" input")));
      await session.step(49, "Then model preview should be ready", () => shouldBe(page, el("model preview"), "ready"));
      await session.step(50, "And \"Impute missing\" input should be hidden", () => shouldBe(page, el("\"Impute missing\" input"), "hidden"));
      await session.step(51, "And \"Predict probability\" input should be visible", () => shouldBe(page, el("\"Predict probability\" input"), "visible"));
      await session.step(52, "And \"Model Engine\" input should be visible", () => shouldBe(page, el("\"Model Engine\" input"), "visible"));
      await session.step(53, "And \"Accuracy\" table row should be visible", () => shouldBe(page, el("\"Accuracy\" table row"), "visible"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A change to a numeric target clears the options ticked, and Predict probability goes", async () => {
      await session.step(57, "When user selects \"AGE\" in Predict input", () => selectIn(page, "AGE", el("Predict input")));
      await session.step(58, "Then editor of Predict input should have text \"AGE\"", () => shouldHaveText(page, el("editor of Predict input"), "AGE"));
      await session.step(59, "And \"Ignore missing\" input should be unchecked", () => shouldBe(page, el("\"Ignore missing\" input"), "unchecked"));
      await session.step(60, "And model preview should be invalid", () => shouldBe(page, el("model preview"), "invalid"));
      await session.step(61, "And \"Model Engine\" input should be absent", () => shouldBe(page, el("\"Model Engine\" input"), "absent"));
      await session.step(62, "When user checks \"Ignore missing\" input", () => check(page, el("\"Ignore missing\" input")));
      await session.step(63, "Then model preview should be ready", () => shouldBe(page, el("model preview"), "ready"));
      await session.step(64, "And \"R squared\" table row should be visible", () => shouldBe(page, el("\"R squared\" table row"), "visible"));
      await session.step(65, "And \"Predict probability\" input should be hidden", () => shouldBe(page, el("\"Predict probability\" input"), "hidden"));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Predict probability trains a regression on the two classes, with a cutoff and a ROC curve", async () => {
      await session.step(69, "Given user opens a table \"readings\" with:", () => openTableOf(page, "readings", [["f1","f2","target"],["1","7","A"],["2","3","A"],["3","9","B"],["4","1","A"],["5","8","B"],["6","2","A"],["7","10","B"],["8","4","B"],["9","6","A"],["10","5","B"]]), [["f1","f2","target"],["1","7","A"],["2","3","A"],["3","9","B"],["4","1","A"],["5","8","B"],["6","2","A"],["7","10","B"],["8","4","B"],["9","6","A"],["10","5","B"]]);
      await session.step(81, "When user picks \"ML > Models > Train Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Train Model..."));
      await session.step(82, "Then the \"Predictive model\" view should be current", () => viewIsCurrent(page, "Predictive model"));
      await session.step(83, "When user selects \"target\" in Predict input", () => selectIn(page, "target", el("Predict input")));
      await session.step(84, "And user clicks on editor of Features input", () => clickOn(page, el("editor of Features input")));
      await session.step(85, "And user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(86, "Then the \"text of cell 1 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"f1\"", () => readingReads(page, "text of cell 1 of __name", el("grid viewer in \"Select columns...\" dialog"), "f1"));
      await session.step(87, "And the \"text of cell 2 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"f2\"", () => readingReads(page, "text of cell 2 of __name", el("grid viewer in \"Select columns...\" dialog"), "f2"));
      await session.step(88, "When user clicks on the \"cell 1 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 1 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(89, "And user clicks on the \"cell 2 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 2 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(90, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(91, "Then model preview should be ready", () => shouldBe(page, el("model preview"), "ready"));
      await session.step(92, "And \"Predict probability\" input should be visible", () => shouldBe(page, el("\"Predict probability\" input"), "visible"));
      await session.step(93, "And \"Positive class cutoff\" input should be absent", () => shouldBe(page, el("\"Positive class cutoff\" input"), "absent"));
      await session.step(94, "And model preview should not contain text \"ROC Curve\"", () => shouldNotContainText(page, el("model preview"), "ROC Curve"));
      await session.step(95, "When user checks \"Predict probability\" input", () => check(page, el("\"Predict probability\" input")));
      await session.step(96, "Then model preview should be ready", () => shouldBe(page, el("model preview"), "ready"));
      await session.step(97, "And \"Positive class cutoff\" input should be visible", () => shouldBe(page, el("\"Positive class cutoff\" input"), "visible"));
      await session.step(98, "And model preview should contain text \"ROC Curve\"", () => shouldContainText(page, el("model preview"), "ROC Curve"));
      await session.step(99, "And \"Eda: Linear Regression\" heading should be visible", () => shouldBe(page, el("\"Eda: Linear Regression\" heading"), "visible"));
      await session.step(100, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(101, "And user enters \"BDD-Probability-{run}\" into Name input in dialog", () => enterInto(page, session.text("BDD-Probability-{run}"), el("Name input in dialog")));
      await session.step(102, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
      await session.step(103, "Then dialog should be absent", () => shouldBe(page, el("dialog"), "absent"));
      await session.step(104, "And 1 predictive model named \"BDD-Probability-{run}\" should be on the server", () => modelsOnServer(page, 1, session.text("BDD-Probability-{run}")));
      await session.step(105, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The probability model applied to new rows writes the class names", async () => {
      await session.step(108, "Given user opens a table \"new readings\" with:", () => openTableOf(page, "new readings", [["f1","f2"],["1","7"],["3","9"],["4","1"],["7","10"],["8","4"],["9","6"]]), [["f1","f2"],["1","7"],["3","9"],["4","1"],["7","10"],["8","4"],["9","6"]]);
      await session.step(116, "When user picks \"ML > Models > Apply Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Apply Model..."));
      await session.step(117, "Then \"Apply predictive model\" dialog should be visible", () => shouldBe(page, el("\"Apply predictive model\" dialog"), "visible"));
      await session.step(118, "When user selects the predictive model \"BDD-Probability-{run}\" in Model input in \"Apply predictive model\" dialog", () => selectModel(page, session.text("BDD-Probability-{run}"), el("Model input in \"Apply predictive model\" dialog")));
      await session.step(119, "Then Inputs input in \"Apply predictive model\" dialog should contain text \"(2/2)\"", () => shouldContainText(page, el("Inputs input in \"Apply predictive model\" dialog"), "(2/2)"));
      await session.step(120, "When user clicks on OK button in \"Apply predictive model\" dialog", () => clickOn(page, el("OK button in \"Apply predictive model\" dialog")));
      await session.step(121, "Then the \"Apply predictive model\" dialog should close", () => dialogCloses(page, "Apply predictive model"));
      await session.step(122, "And a new column \"target\" should have been added", () => newColumnNamed(page, "target"));
      await session.step(123, "And \"target\" column should have no missing values", () => columnComplete(page, "target"));
      await session.step(125, "And the \"target\" cell of row 1 should be displayed as \"A\"", () => displayedInRow(page, "target", 1, "A"));
      await session.step(126, "And the \"target\" cell of row 4 should be displayed as \"B\"", () => displayedInRow(page, "target", 4, "B"));
      await session.step(127, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("One-hot encoding trains on two Yes/No features, and the model applies to a fresh table", async () => {
      await session.step(130, "Given user opens a table \"answers\" with:", () => openTableOf(page, "answers", [["featureA","featureB","target"],["Yes","Yes","5.00"],["No","Yes","3.07"],["Yes","No","2.14"],["No","No","0.21"],["Yes","Yes","5.28"],["No","Yes","3.35"],["Yes","No","2.42"],["No","No","0.49"],["Yes","Yes","5.56"],["No","Yes","3.63"],["Yes","No","2.70"],["No","No","0.77"]]), [["featureA","featureB","target"],["Yes","Yes","5.00"],["No","Yes","3.07"],["Yes","No","2.14"],["No","No","0.21"],["Yes","Yes","5.28"],["No","Yes","3.35"],["Yes","No","2.42"],["No","No","0.49"],["Yes","Yes","5.56"],["No","Yes","3.63"],["Yes","No","2.70"],["No","No","0.77"]]);
      await session.step(144, "When user picks \"ML > Models > Train Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Train Model..."));
      await session.step(145, "Then the \"Predictive model\" view should be current", () => viewIsCurrent(page, "Predictive model"));
      await session.step(146, "When user selects \"target\" in Predict input", () => selectIn(page, "target", el("Predict input")));
      await session.step(147, "And user clicks on editor of Features input", () => clickOn(page, el("editor of Features input")));
      await session.step(148, "And user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(149, "Then the \"text of cell 1 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"featureA\"", () => readingReads(page, "text of cell 1 of __name", el("grid viewer in \"Select columns...\" dialog"), "featureA"));
      await session.step(150, "And the \"text of cell 2 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"featureB\"", () => readingReads(page, "text of cell 2 of __name", el("grid viewer in \"Select columns...\" dialog"), "featureB"));
      await session.step(151, "When user clicks on the \"cell 1 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 1 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(152, "And user clicks on the \"cell 2 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 2 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(153, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(154, "Then model preview should contain text \"Columns 'featureA, featureB' are categorical.\"", () => shouldContainText(page, el("model preview"), "Columns 'featureA, featureB' are categorical."));
      await session.step(155, "And model preview should be invalid", () => shouldBe(page, el("model preview"), "invalid"));
      await session.step(156, "And \"Model Engine\" input should be absent", () => shouldBe(page, el("\"Model Engine\" input"), "absent"));
      await session.step(157, "When user checks \"One-hot encoding\" input", () => check(page, el("\"One-hot encoding\" input")));
      await session.step(158, "Then model preview should be ready", () => shouldBe(page, el("model preview"), "ready"));
      await session.step(159, "And \"Model Engine\" input should be visible", () => shouldBe(page, el("\"Model Engine\" input"), "visible"));
      await session.step(160, "When user selects \"Eda: Linear Regression\" in \"Model Engine\" input", () => selectIn(page, "Eda: Linear Regression", el("\"Model Engine\" input")));
      await session.step(161, "Then model preview should be ready", () => shouldBe(page, el("model preview"), "ready"));
      await session.step(162, "And \"Eda: Linear Regression\" heading should be visible", () => shouldBe(page, el("\"Eda: Linear Regression\" heading"), "visible"));
      await session.step(163, "And \"R squared\" table row should be visible", () => shouldBe(page, el("\"R squared\" table row"), "visible"));
      await session.step(164, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(165, "And user enters \"BDD-OneHot-{run}\" into Name input in dialog", () => enterInto(page, session.text("BDD-OneHot-{run}"), el("Name input in dialog")));
      await session.step(166, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
      await session.step(167, "Then dialog should be absent", () => shouldBe(page, el("dialog"), "absent"));
      await session.step(168, "And 1 predictive model named \"BDD-OneHot-{run}\" should be on the server", () => modelsOnServer(page, 1, session.text("BDD-OneHot-{run}")));
      await session.step(169, "Given user opens a table \"new answers\" with:", () => openTableOf(page, "new answers", [["featureA","featureB"],["No","No"],["Yes","No"],["No","Yes"],["Yes","Yes"]]), [["featureA","featureB"],["No","No"],["Yes","No"],["No","Yes"],["Yes","Yes"]]);
      await session.step(175, "When user picks \"ML > Models > Apply Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Apply Model..."));
      await session.step(176, "And user selects the predictive model \"BDD-OneHot-{run}\" in Model input in \"Apply predictive model\" dialog", () => selectModel(page, session.text("BDD-OneHot-{run}"), el("Model input in \"Apply predictive model\" dialog")));
      await session.step(177, "Then Inputs input in \"Apply predictive model\" dialog should contain text \"(2/2)\"", () => shouldContainText(page, el("Inputs input in \"Apply predictive model\" dialog"), "(2/2)"));
      await session.step(178, "When user clicks on OK button in \"Apply predictive model\" dialog", () => clickOn(page, el("OK button in \"Apply predictive model\" dialog")));
      await session.step(179, "Then the \"Apply predictive model\" dialog should close", () => dialogCloses(page, "Apply predictive model"));
      await session.step(180, "And a new column \"target\" should have been added", () => newColumnNamed(page, "target"));
      await session.step(181, "And \"target\" column should have no missing values", () => columnComplete(page, "target"));
      await session.step(183, "And the \"target\" cell of row 2 should be displayed as \"2.42\"", () => displayedInRow(page, "target", 2, "2.42"));
      await session.step(184, "And the \"target\" cell of row 3 should be displayed as \"3.35\"", () => displayedInRow(page, "target", 3, "3.35"));
      await session.step(185, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(186, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
