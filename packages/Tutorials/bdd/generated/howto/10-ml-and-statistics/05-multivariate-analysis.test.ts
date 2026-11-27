/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/10-ml-and-statistics/05-multivariate-analysis.feature
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
import {clickOn, enterInto, hoverText, selectIn, shouldBe, shouldContainText, typeInto, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · ML", () => {
  const session = feature(test, "features/howto/10-ml-and-statistics/05-multivariate-analysis.feature", import.meta.url);
  test("Find what drives a value with multivariate analysis", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(19, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(22, "When the video says \"A car's price, and fifteen specs that might explain it\" at the bottom", () => videoSaysAt(page, "A car's price, and fifteen specs that might explain it", "bottom"));
    await session.step(24, "And user picks \"ML > Analyze > Multivariate Analysis...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > Multivariate Analysis..."), undefined, "Multivariate Analysis is under ML > Analyze");
    await session.step(25, "Then \"Multivariate Analysis (PLS)\" dialog should be visible", () => shouldBe(page, el("\"Multivariate Analysis (PLS)\" dialog"), "visible"));
    await session.step(27, "When user selects \"price\" in Predict input in \"Multivariate Analysis (PLS)\" dialog", () => selectIn(page, "price", el("Predict input in \"Multivariate Analysis (PLS)\" dialog")), undefined, "Predict: the value to explain");
    await session.step(29, "And user clicks on editor of Using input in \"Multivariate Analysis (PLS)\" dialog", () => clickOn(page, el("editor of Using input in \"Multivariate Analysis (PLS)\" dialog")), undefined, "Using: the columns that might explain it");
    await session.step(30, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(31, "When user clicks on \"All\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"All\" link in \"Select columns...\" dialog")));
    await session.step(33, "And user types \"price\" into Search input in \"Select columns...\" dialog", () => typeInto(page, "price", el("Search input in \"Select columns...\" dialog")), undefined, "Leave out price itself: it is the value to explain");
    await session.step(34, "And user toggles the \"price\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "price", el("\"Select columns...\" dialog")));
    await session.step(35, "Then \"Select columns...\" dialog should contain text \"15 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "15 checked"));
    await session.step(36, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(38, "And user enters \"3\" into Components input in \"Multivariate Analysis (PLS)\" dialog", () => enterInto(page, "3", el("Components input in \"Multivariate Analysis (PLS)\" dialog")), undefined, "Components: how many hidden factors to fit");
    await session.step(40, "And user selects \"model\" in Names input in \"Multivariate Analysis (PLS)\" dialog", () => selectIn(page, "model", el("Names input in \"Multivariate Analysis (PLS)\" dialog")), undefined, "Names: label the points with the car model");
    await session.step(41, "And user clicks on RUN button in \"Multivariate Analysis (PLS)\" dialog", () => clickOn(page, el("RUN button in \"Multivariate Analysis (PLS)\" dialog")));
    await session.step(42, "Then the table should have a column \"price (predicted)\"", () => hasColumn(page, "price (predicted)"));
    await session.step(43, "And the open tableview should have 3 scatter plot viewers", () => viewerCount(page, 3, "scatter plot"));
    await session.step(44, "And the open tableview should have 3 bar chart viewers", () => viewerCount(page, 3, "bar chart"));
    await session.step(46, "When user hovers over the text \"Observed vs. Predicted\" in \"Observed vs. Predicted\" dock panel", () => hoverText(page, "Observed vs. Predicted", el("\"Observed vs. Predicted\" dock panel")), undefined, "Observed vs. Predicted: close to the line means a good fit");
    await session.step(48, "And user hovers over the text \"Variable Importance\" in \"Variable Importance\" dock panel", () => hoverText(page, "Variable Importance", el("\"Variable Importance\" dock panel")), undefined, "Variable Importance: the specs that matter most for price");
    await session.step(49, "When the video says \"Loadings: specs close together move together\" at the bottom", () => videoSaysAt(page, "Loadings: specs close together move together", "bottom"));
    await session.step(50, "And user hovers over the text \"Loadings\" in \"Loadings\" dock panel", () => hoverText(page, "Loadings", el("\"Loadings\" dock panel")));
    await session.step(51, "When the video says \"Scores: each car placed by the fitted factors\" at the bottom", () => videoSaysAt(page, "Scores: each car placed by the fitted factors", "bottom"));
    await session.step(52, "And user hovers over the text \"Scores\" in \"Scores\" dock panel", () => hoverText(page, "Scores", el("\"Scores\" dock panel")));
    await session.step(54, "And user hovers over the text \"Explained Variance\" in \"Explained Variance\" dock panel", () => hoverText(page, "Explained Variance", el("\"Explained Variance\" dock panel")), undefined, "Explained Variance: share of price explained by 1, 2 and 3 components");
    await session.step(55, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(56, "And no errors should have been logged", () => noErrors(page));
  });
});
