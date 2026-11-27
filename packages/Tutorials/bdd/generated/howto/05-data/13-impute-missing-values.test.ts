/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/05-data/13-impute-missing-values.feature
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
import {clickOn, shouldBe, shouldContainText, uncheck, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnIncomplete, hasColumn, makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {columnListExactly, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Data", () => {
  const session = feature(test, "features/howto/05-data/13-impute-missing-values.feature", import.meta.url);
  test("Fill in missing values", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(19, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(22, "And user makes row 2222 current", () => makeRowCurrent(page, 2222));
    await session.step(23, "Then \"HEIGHT\" column should have missing values", () => columnIncomplete(page, "HEIGHT"));
    await session.step(24, "And \"WEIGHT\" column should have missing values", () => columnIncomplete(page, "WEIGHT"));
    await session.step(25, "And \"AGE\" column should have missing values", () => columnIncomplete(page, "AGE"));
    await session.step(26, "When the video says \"Empty cells: 751 heights, one weight and one age are missing\"", () => videoSays(page, "Empty cells: 751 heights, one weight and one age are missing"));
    await session.step(28, "And user picks \"ML > Impute Missing Values...\" from the top menu", () => pickFromTopMenu(page, "ML > Impute Missing Values..."), undefined, "Impute Missing Values fills gaps from similar rows");
    await session.step(29, "Then \"k-NN Imputation\" dialog should be visible", () => shouldBe(page, el("\"k-NN Imputation\" dialog"), "visible"));
    await session.step(31, "When user clicks on Impute input in \"k-NN Imputation\" dialog", () => clickOn(page, el("Impute input in \"k-NN Imputation\" dialog")), undefined, "Impute: the columns with gaps, all three ticked");
    await session.step(32, "Then the column list of \"Select columns...\" dialog should be exactly \"AGE, HEIGHT, WEIGHT\"", () => columnListExactly(page, el("\"Select columns...\" dialog"), "AGE, HEIGHT, WEIGHT"));
    await session.step(33, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(35, "And user clicks on Using input in \"k-NN Imputation\" dialog", () => clickOn(page, el("Using input in \"k-NN Imputation\" dialog")), undefined, "Using: the columns that decide which rows are alike");
    await session.step(37, "And user toggles the \"USUBJID\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "USUBJID", el("\"Select columns...\" dialog")), undefined, "An ID says nothing about a row: leave it out");
    await session.step(38, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(39, "Then Using input in \"k-NN Imputation\" dialog should contain text \"(9)\"", () => shouldContainText(page, el("Using input in \"k-NN Imputation\" dialog"), "(9)"));
    await session.step(41, "When user unchecks \"In-place\" input in \"k-NN Imputation\" dialog", () => uncheck(page, el("\"In-place\" input in \"k-NN Imputation\" dialog")), undefined, "Untick In-place to keep the originals next to filled copies");
    await session.step(43, "And user clicks on RUN button in \"k-NN Imputation\" dialog", () => clickOn(page, el("RUN button in \"k-NN Imputation\" dialog")), undefined, "RUN fills each gap from its 4 nearest rows");
    await session.step(44, "Then \"k-NN Imputation\" dialog should be hidden", () => shouldBe(page, el("\"k-NN Imputation\" dialog"), "hidden"));
    await session.step(45, "And the table should have a column \"HEIGHT(copy)\"", () => hasColumn(page, "HEIGHT(copy)"));
    await session.step(46, "And \"HEIGHT(copy)\" column should have no missing values", () => columnComplete(page, "HEIGHT(copy)"));
    await session.step(47, "And \"WEIGHT(copy)\" column should have no missing values", () => columnComplete(page, "WEIGHT(copy)"));
    await session.step(48, "And \"AGE(copy)\" column should have no missing values", () => columnComplete(page, "AGE(copy)"));
    await session.step(49, "And \"HEIGHT\" column should have missing values", () => columnIncomplete(page, "HEIGHT"));
    await session.step(50, "When the video says \"The copies have no gaps; the originals are untouched\" at the bottom", () => videoSaysAt(page, "The copies have no gaps; the originals are untouched", "bottom"));
    await session.step(51, "And user hovers over the \"cell 2221 of HEIGHT(copy)\" area of grid", () => hoverArea(page, "cell 2221 of HEIGHT(copy)", el("grid")));
    await session.step(52, "And user hovers over the \"cell 2222 of WEIGHT(copy)\" area of grid", () => hoverArea(page, "cell 2222 of WEIGHT(copy)", el("grid")));
    await session.step(53, "Then no errors should have been logged", () => noErrors(page));
  });
});
