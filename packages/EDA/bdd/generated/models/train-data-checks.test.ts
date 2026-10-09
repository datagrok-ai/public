/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/models/train-data-checks.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [ml.menu.models.train-model]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, selectIn, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {closeCurrentView, openTableOf, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Train Model view checks the data before it trains", () => {
  const session = feature(test, "features/models/train-data-checks.feature", import.meta.url);
  test("The Train Model view checks the data before it trains", {tag: ["@journey", "@eda", "@realizes:ml.menu.models.train-model"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await run.scenario("Class imbalance is reported and the model still trains", async () => {
      await session.step(33, "Given user opens a table \"imbalanced\" with:", () => openTableOf(page, "imbalanced", [["f1","f2","target"],["1","7","A"],["2","3","A"],["3","9","A"],["4","1","A"],["5","8","A"],["6","2","A"],["7","10","A"],["8","4","A"],["9","6","B"],["10","5","B"]]), [["f1","f2","target"],["1","7","A"],["2","3","A"],["3","9","A"],["4","1","A"],["5","8","A"],["6","2","A"],["7","10","A"],["8","4","A"],["9","6","B"],["10","5","B"]]);
      await session.step(45, "When user picks \"ML > Models > Train Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Train Model..."));
      await session.step(46, "Then the \"Predictive model\" view should be current", () => viewIsCurrent(page, "Predictive model"));
      await session.step(47, "When user selects \"target\" in Predict input", () => selectIn(page, "target", el("Predict input")));
      await session.step(48, "And user clicks on editor of Features input", () => clickOn(page, el("editor of Features input")));
      await session.step(49, "And user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(50, "And user clicks on the \"cell 1 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 1 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(51, "And user clicks on the \"cell 2 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 2 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(52, "Then \"Select columns...\" dialog should contain text \"2 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "2 checked"));
      await session.step(53, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(54, "Then editor of Features input should contain text \"(2)\"", () => shouldContainText(page, el("editor of Features input"), "(2)"));
      await session.step(55, "And model preview should be ready", () => shouldBe(page, el("model preview"), "ready"));
      await session.step(56, "And model preview should contain text \"Some columns contain class imbalance\"", () => shouldContainText(page, el("model preview"), "Some columns contain class imbalance"));
      await session.step(57, "And model preview should contain text \"A (1.600)\"", () => shouldContainText(page, el("model preview"), "A (1.600)"));
      await session.step(58, "And model preview should contain text \"B (0.400)\"", () => shouldContainText(page, el("model preview"), "B (0.400)"));
      await session.step(59, "And \"Model Engine\" input should be visible", () => shouldBe(page, el("\"Model Engine\" input"), "visible"));
      await session.step(60, "And Save button should be enabled", () => shouldBe(page, el("Save button"), "enabled"));
      await session.step(61, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
      await session.step(63, "When user closes the current view", () => closeCurrentView(page));
    });
    await run.scenario("A categorical feature is reported, and One-hot encoding lets the model train", async () => {
      await session.step(66, "Given user opens a table \"colors\" with:", () => openTableOf(page, "colors", [["cat","num","target"],["red","1","3.1"],["green","2","1.2"],["blue","3","4.4"],["red","4","2.5"],["green","5","5.3"],["blue","6","1.7"],["red","7","3.9"],["green","8","2.2"],["blue","9","4.8"]]), [["cat","num","target"],["red","1","3.1"],["green","2","1.2"],["blue","3","4.4"],["red","4","2.5"],["green","5","5.3"],["blue","6","1.7"],["red","7","3.9"],["green","8","2.2"],["blue","9","4.8"]]);
      await session.step(77, "When user picks \"ML > Models > Train Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Train Model..."));
      await session.step(78, "And user selects \"target\" in Predict input", () => selectIn(page, "target", el("Predict input")));
      await session.step(79, "And user clicks on editor of Features input", () => clickOn(page, el("editor of Features input")));
      await session.step(80, "And user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(81, "And user clicks on the \"cell 1 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 1 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(82, "And user clicks on the \"cell 2 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 2 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(83, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(84, "Then editor of Features input should contain text \"(2)\"", () => shouldContainText(page, el("editor of Features input"), "(2)"));
      await session.step(85, "And model preview should contain text \"Column 'cat' is categorical. Most models require converting them to numerical.\"", () => shouldContainText(page, el("model preview"), "Column 'cat' is categorical. Most models require converting them to numerical."));
      await session.step(86, "And model preview should be invalid", () => shouldBe(page, el("model preview"), "invalid"));
      await session.step(87, "And \"One-hot encoding\" input should be visible", () => shouldBe(page, el("\"One-hot encoding\" input"), "visible"));
      await session.step(88, "And \"Model Engine\" input should be absent", () => shouldBe(page, el("\"Model Engine\" input"), "absent"));
      await session.step(89, "And Save button should be disabled", () => shouldBe(page, el("Save button"), "disabled"));
      await session.step(90, "When user checks \"One-hot encoding\" input", () => check(page, el("\"One-hot encoding\" input")));
      await session.step(91, "Then model preview should be ready", () => shouldBe(page, el("model preview"), "ready"));
      await session.step(92, "And \"Model Engine\" input should be visible", () => shouldBe(page, el("\"Model Engine\" input"), "visible"));
      await session.step(93, "And model preview should not contain text \"is categorical\"", () => shouldNotContainText(page, el("model preview"), "is categorical"));
      await session.step(94, "And \"R squared\" table row should be visible", () => shouldBe(page, el("\"R squared\" table row"), "visible"));
      await session.step(95, "And Save button should be enabled", () => shouldBe(page, el("Save button"), "enabled"));
      await session.step(96, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
      await session.step(98, "When user closes the current view", () => closeCurrentView(page));
    });
    await run.scenario("An identifier-like column is reported, and the view offers to skip it", async () => {
      await session.step(101, "Given user opens a table \"identifiers\" with:", () => openTableOf(page, "identifiers", [["id_like","feature1","target"],["row_1","1","3.1"],["row_2","2","1.2"],["row_3","3","4.4"],["row_4","4","2.5"],["row_5","5","5.3"],["row_6","6","1.7"]]), [["id_like","feature1","target"],["row_1","1","3.1"],["row_2","2","1.2"],["row_3","3","4.4"],["row_4","4","2.5"],["row_5","5","5.3"],["row_6","6","1.7"]]);
      await session.step(109, "When user picks \"ML > Models > Train Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Train Model..."));
      await session.step(110, "And user selects \"target\" in Predict input", () => selectIn(page, "target", el("Predict input")));
      await session.step(111, "And user clicks on editor of Features input", () => clickOn(page, el("editor of Features input")));
      await session.step(112, "And user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(113, "And user clicks on the \"cell 1 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 1 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(114, "And user clicks on the \"cell 2 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 2 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(115, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(116, "Then editor of Features input should contain text \"(2)\"", () => shouldContainText(page, el("editor of Features input"), "(2)"));
      await session.step(117, "And model preview should contain text \"Column 'id_like' contains\"", () => shouldContainText(page, el("model preview"), "Column 'id_like' contains"));
      await session.step(118, "And model preview should contain text \"too many unique categories.\"", () => shouldContainText(page, el("model preview"), "too many unique categories."));
      await session.step(119, "And model preview should be invalid", () => shouldBe(page, el("model preview"), "invalid"));
      await session.step(120, "And \"Skip unique categories\" input should be visible", () => shouldBe(page, el("\"Skip unique categories\" input"), "visible"));
      await session.step(121, "And \"Model Engine\" input should be absent", () => shouldBe(page, el("\"Model Engine\" input"), "absent"));
      await session.step(122, "And Save button should be disabled", () => shouldBe(page, el("Save button"), "disabled"));
      await session.step(123, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(124, "And no errors should have been logged", () => noErrors(page));
      await session.step(125, "When user closes the current view", () => closeCurrentView(page));
    });
    await run.scenario("Highly correlated features are named in pairs and the model still trains", async () => {
      await session.step(128, "Given user opens a table \"correlated\" with:", () => openTableOf(page, "correlated", [["feat_a","feat_b","target"],["1","1.1","3.1"],["2","2.3","1.2"],["3","2.9","4.4"],["4","4.2","2.5"],["5","5.1","5.3"],["6","6.4","1.7"],["7","6.8","3.9"],["8","8.3","2.2"]]), [["feat_a","feat_b","target"],["1","1.1","3.1"],["2","2.3","1.2"],["3","2.9","4.4"],["4","4.2","2.5"],["5","5.1","5.3"],["6","6.4","1.7"],["7","6.8","3.9"],["8","8.3","2.2"]]);
      await session.step(138, "When user picks \"ML > Models > Train Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Train Model..."));
      await session.step(139, "And user selects \"target\" in Predict input", () => selectIn(page, "target", el("Predict input")));
      await session.step(140, "And user clicks on editor of Features input", () => clickOn(page, el("editor of Features input")));
      await session.step(141, "And user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(142, "And user clicks on the \"cell 1 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 1 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(143, "And user clicks on the \"cell 2 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 2 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(144, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(145, "Then editor of Features input should contain text \"(2)\"", () => shouldContainText(page, el("editor of Features input"), "(2)"));
      await session.step(146, "And model preview should be ready", () => shouldBe(page, el("model preview"), "ready"));
      await session.step(147, "And model preview should contain text \"Columns are highly correlated\"", () => shouldContainText(page, el("model preview"), "Columns are highly correlated"));
      await session.step(148, "And model preview should contain text \"feat_a <> feat_b\"", () => shouldContainText(page, el("model preview"), "feat_a <> feat_b"));
      await session.step(149, "And Save button should be enabled", () => shouldBe(page, el("Save button"), "enabled"));
      await session.step(150, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(151, "And no errors should have been logged", () => noErrors(page));
      await session.step(152, "When user closes the current view", () => closeCurrentView(page));
    });
    await run.scenario("Missing values in a numeric target block training until Ignore missing is ticked", async () => {
      await session.step(156, "Given user opens a table \"gaps\" with:", () => openTableOf(page, "gaps", [["X","Y"],["1","2"],["2","4"],["3",""],["4","8"],["5","10"],["6",""],["7","14"],["8","16"],["9","18"],["10","20"]]), [["X","Y"],["1","2"],["2","4"],["3",""],["4","8"],["5","10"],["6",""],["7","14"],["8","16"],["9","18"],["10","20"]]);
      await session.step(168, "When user picks \"ML > Models > Train Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Train Model..."));
      await session.step(169, "And user selects \"Y\" in Predict input", () => selectIn(page, "Y", el("Predict input")));
      await session.step(170, "And user clicks on editor of Features input", () => clickOn(page, el("editor of Features input")));
      await session.step(171, "And user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(172, "And user clicks on the \"cell 1 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 1 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(173, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(174, "Then editor of Features input should contain text \"(1)\"", () => shouldContainText(page, el("editor of Features input"), "(1)"));
      await session.step(175, "And model preview should contain text \"Column 'Y' contains missing values.\"", () => shouldContainText(page, el("model preview"), "Column 'Y' contains missing values."));
      await session.step(176, "And model preview should be invalid", () => shouldBe(page, el("model preview"), "invalid"));
      await session.step(177, "And \"Ignore missing\" input should be visible", () => shouldBe(page, el("\"Ignore missing\" input"), "visible"));
      await session.step(178, "And \"Impute missing\" input should be visible", () => shouldBe(page, el("\"Impute missing\" input"), "visible"));
      await session.step(179, "And \"Model Engine\" input should be absent", () => shouldBe(page, el("\"Model Engine\" input"), "absent"));
      await session.step(180, "And Save button should be disabled", () => shouldBe(page, el("Save button"), "disabled"));
      await session.step(181, "When user checks \"Ignore missing\" input", () => check(page, el("\"Ignore missing\" input")));
      await session.step(182, "Then model preview should be ready", () => shouldBe(page, el("model preview"), "ready"));
      await session.step(183, "And \"Impute missing\" input should be hidden", () => shouldBe(page, el("\"Impute missing\" input"), "hidden"));
      await session.step(184, "And \"Model Engine\" input should be visible", () => shouldBe(page, el("\"Model Engine\" input"), "visible"));
      await session.step(185, "And model preview should not contain text \"contains missing values\"", () => shouldNotContainText(page, el("model preview"), "contains missing values"));
      await session.step(186, "And \"R squared\" table row should be visible", () => shouldBe(page, el("\"R squared\" table row"), "visible"));
      await session.step(187, "And Save button should be enabled", () => shouldBe(page, el("Save button"), "enabled"));
      await session.step(188, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(189, "And no errors should have been logged", () => noErrors(page));
      await session.step(190, "When user closes the current view", () => closeCurrentView(page));
    });
    await run.scenario("Missing labels in a categorical target block training too", async () => {
      await session.step(194, "Given user opens a table \"unlabelled\" with:", () => openTableOf(page, "unlabelled", [["f1","f2","label"],["1","7","A"],["2","3","B"],["3","9",""],["4","1","A"],["5","8","B"],["6","2","A"],["7","10","B"],["8","4",""]]), [["f1","f2","label"],["1","7","A"],["2","3","B"],["3","9",""],["4","1","A"],["5","8","B"],["6","2","A"],["7","10","B"],["8","4",""]]);
      await session.step(204, "When user picks \"ML > Models > Train Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Train Model..."));
      await session.step(205, "And user selects \"label\" in Predict input", () => selectIn(page, "label", el("Predict input")));
      await session.step(206, "And user clicks on editor of Features input", () => clickOn(page, el("editor of Features input")));
      await session.step(207, "And user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(208, "And user clicks on the \"cell 1 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 1 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(209, "And user clicks on the \"cell 2 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 2 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(210, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(211, "Then editor of Features input should contain text \"(2)\"", () => shouldContainText(page, el("editor of Features input"), "(2)"));
      await session.step(212, "And model preview should contain text \"Column 'label' contains missing values.\"", () => shouldContainText(page, el("model preview"), "Column 'label' contains missing values."));
      await session.step(213, "And model preview should be invalid", () => shouldBe(page, el("model preview"), "invalid"));
      await session.step(214, "And \"Model Engine\" input should be absent", () => shouldBe(page, el("\"Model Engine\" input"), "absent"));
      await session.step(215, "And Save button should be disabled", () => shouldBe(page, el("Save button"), "disabled"));
      await session.step(216, "When user checks \"Ignore missing\" input", () => check(page, el("\"Ignore missing\" input")));
      await session.step(217, "Then model preview should be ready", () => shouldBe(page, el("model preview"), "ready"));
      await session.step(218, "And \"Model Engine\" input should be visible", () => shouldBe(page, el("\"Model Engine\" input"), "visible"));
      await session.step(219, "And model preview should not contain text \"contains missing values\"", () => shouldNotContainText(page, el("model preview"), "contains missing values"));
      await session.step(220, "And Save button should be enabled", () => shouldBe(page, el("Save button"), "enabled"));
      await session.step(221, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(222, "And no errors should have been logged", () => noErrors(page));
      await session.step(223, "When user closes the current view", () => closeCurrentView(page));
    });
    run.finish();
  });
});
