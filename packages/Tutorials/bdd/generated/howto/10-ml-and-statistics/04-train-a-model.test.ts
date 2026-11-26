/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/10-ml-and-statistics/04-train-a-model.feature
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
import {forgeSweep, newTrainingShown, rememberShownTraining, trainingShown} from '../../../bindings/howto-ml.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, expand, hoverOver, hoverText, selectIn, shouldBe, shouldContainText, shouldHaveValue, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, distinctValues, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelClosed, clickViewTab, contextPanelClosed, currentViewShows, dialogCloses, openDatasetRowsAs, openTableOf, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · ML", () => {
  const session = feature(test, "features/howto/10-ml-and-statistics/04-train-a-model.feature", import.meta.url);
  test("Train a model and apply it", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And simple mode is off", () => simpleModeOff(page));
    await session.step(26, "And the Forge model \"Iris species\" and the training runs of the table \"iris flowers\" are removed now and when the feature ends", () => forgeSweep(page, "Iris species", "iris flowers"));
    await session.step(27, "And user opens a table \"new flowers\" with:", () => openTableOf(page, "new flowers", [["Sepal.Length","Sepal.Width","Petal.Length","Petal.Width"],["5.1","3.5","1.4","0.2"],["7.0","3.2","4.7","1.4"],["6.3","3.3","6.0","2.5"],["5.4","3.7","1.5","0.2"],["5.0","2.0","3.5","1.0"],["6.5","3.2","5.1","2.0"],["5.4","3.4","1.7","0.2"],["5.5","2.4","3.8","1.1"],["6.9","3.2","5.7","2.3"]]), [["Sepal.Length","Sepal.Width","Petal.Length","Petal.Width"],["5.1","3.5","1.4","0.2"],["7.0","3.2","4.7","1.4"],["6.3","3.3","6.0","2.5"],["5.4","3.7","1.5","0.2"],["5.0","2.0","3.5","1.0"],["6.5","3.2","5.1","2.0"],["5.4","3.4","1.7","0.2"],["5.5","2.4","3.8","1.1"],["6.9","3.2","5.7","2.3"]]);
    await session.step(38, "And user opens iris dataset keeping the first 150 rows as \"iris flowers\"", () => openDatasetRowsAs(page, ds("iris"), 150, "iris flowers"));
    await session.step(39, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(40, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(41, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(42, "When the video says \"150 flowers: four measurements and the species\" at the top", () => videoSaysAt(page, "150 flowers: four measurements and the species", "top"));
    await session.step(44, "And user picks \"ML > Models > Train Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Train Model..."), undefined, "Train Model is under ML > Models");
    await session.step(45, "Then the Train view should show a finished training", () => trainingShown(page));
    await session.step(47, "When user hovers over editor of Target input", () => hoverOver(page, el("editor of Target input")), undefined, "Target: the column to predict");
    await session.step(49, "And user hovers over editor of Features input", () => hoverOver(page, el("editor of Features input")), undefined, "Features: the columns it learns from");
    await session.step(51, "And user hovers over editor of Method input", () => hoverOver(page, el("editor of Method input")), undefined, "Method: XGBoost here; training runs in your browser");
    await session.step(53, "And user hovers over title of \"Performance\" results panel", () => hoverOver(page, el("title of \"Performance\" results panel")), undefined, "It trains right away. Performance: scores on held-out rows");
    await session.step(54, "Then \"Performance\" results panel should contain text \"5-fold cross-validation\"", () => shouldContainText(page, el("\"Performance\" results panel"), "5-fold cross-validation"));
    await session.step(56, "When user hovers over title of \"Confusions\" results panel", () => hoverOver(page, el("title of \"Confusions\" results panel")), undefined, "Confusions: predicted against actual, on the rows it learned from");
    await session.step(58, "And user remembers the training the Train view shows", () => rememberShownTraining(page), undefined, "Pick another method and it trains again");
    await session.step(59, "And user selects \"Random Forest\" in Method input", () => selectIn(page, "Random Forest", el("Method input")));
    await session.step(60, "Then the Train view should show a new training", () => newTrainingShown(page));
    await session.step(62, "When user hovers over title of \"Runs\" results panel", () => hoverOver(page, el("title of \"Runs\" results panel")), undefined, "Runs: every training so far, to compare");
    await session.step(64, "When user clicks on Save button", () => clickOn(page, el("Save button")), undefined, "Save keeps the model under a name");
    await session.step(65, "Then \"Save model\" dialog should be visible", () => shouldBe(page, el("\"Save model\" dialog"), "visible"));
    await session.step(66, "When user enters \"Iris species\" into Name input in \"Save model\" dialog", () => enterInto(page, "Iris species", el("Name input in \"Save model\" dialog")));
    await session.step(67, "And user clicks on OK button in \"Save model\" dialog", () => clickOn(page, el("OK button in \"Save model\" dialog")));
    await session.step(68, "Then the \"Save model\" dialog should close", () => dialogCloses(page, "Save model"));
    await session.step(70, "When user clicks on browse tab", () => clickOn(page, el("browse tab")), undefined, "Saved models are in Browse > Platform > Machine Learning");
    await session.step(71, "And user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(72, "And user clicks on \"Machine Learning\" tree node inside browse tree", () => clickOn(page, el("\"Machine Learning\" tree node inside browse tree")));
    await session.step(73, "Then the current view should show the text \"Iris species\"", () => currentViewShows(page, "Iris species"));
    await session.step(75, "When user hovers over the text \"Iris species\" in gallery", () => hoverText(page, "Iris species", el("gallery")), undefined, "Your model's card: what it predicts, from what, and how");
    await session.step(77, "And user clicks on the tab of the \"new flowers\" view", () => clickViewTab(page, "new flowers"), undefined, "Now apply it: a table of nine new flowers");
    await session.step(78, "Then the \"new flowers\" view should be current", () => viewIsCurrent(page, "new flowers"));
    await session.step(79, "When the video says \"Measured, but the species is not known yet\" at the bottom", () => videoSaysAt(page, "Measured, but the species is not known yet", "bottom"));
    await session.step(81, "And user picks \"ML > Models > Apply Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Apply Model..."), undefined, "Apply Model is under ML > Models");
    await session.step(82, "Then \"Apply predictive model\" dialog should be visible", () => shouldBe(page, el("\"Apply predictive model\" dialog"), "visible"));
    await session.step(84, "When user hovers over Model input in \"Apply predictive model\" dialog", () => hoverOver(page, el("Model input in \"Apply predictive model\" dialog")), undefined, "Model: the one just saved is already picked");
    await session.step(85, "Then Model input in \"Apply predictive model\" dialog should have value \"Iris species\"", () => shouldHaveValue(page, el("Model input in \"Apply predictive model\" dialog"), "Iris species"));
    await session.step(87, "When user expands Columns collapsible group in \"Apply predictive model\" dialog", () => expand(page, el("Columns collapsible group in \"Apply predictive model\" dialog")), undefined, "Columns: each input of the model gets the column of the same name");
    await session.step(88, "And user hovers over editor of \"Petal.Length\" input in \"Apply predictive model\" dialog", () => hoverOver(page, el("editor of \"Petal.Length\" input in \"Apply predictive model\" dialog")));
    await session.step(90, "And user clicks on OK button in \"Apply predictive model\" dialog", () => clickOn(page, el("OK button in \"Apply predictive model\" dialog")), undefined, "OK adds the prediction as a new column");
    await session.step(91, "Then \"Apply predictive model\" dialog should be hidden", () => shouldBe(page, el("\"Apply predictive model\" dialog"), "hidden"));
    await session.step(92, "And the table should have a column \"Species (predicted)\"", () => hasColumn(page, "Species (predicted)"));
    await session.step(93, "And \"Species (predicted)\" column should have no missing values", () => columnComplete(page, "Species (predicted)"));
    await session.step(94, "And \"Species (predicted)\" column should have at least 3 distinct values", () => distinctValues(page, "Species (predicted)", 3));
    await session.step(96, "When user hovers over the \"cell 2 of Species (predicted)\" area of grid", () => hoverArea(page, "cell 2 of Species (predicted)", el("grid")), undefined, "Each flower gets the species the model expects from its measurements");
    await session.step(97, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(98, "And no errors should have been logged", () => noErrors(page));
  });
});
