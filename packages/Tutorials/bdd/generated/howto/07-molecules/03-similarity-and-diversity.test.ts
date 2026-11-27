/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/07-molecules/03-similarity-and-diversity.feature
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
import {clickOn, shouldBe, shouldHaveValue, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnMatching, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, openDataset, packageInstalled, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverFirstFreeArea, noErrors, painted, propertyShouldBe, readingNotAsRemembered, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {clickNthCardBesideDrawing, currentRowIsClickedCard, gridShowsColumns, pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Molecules", () => {
  const session = feature(test, "features/howto/07-molecules/03-similarity-and-diversity.feature", import.meta.url);
  test("Find similar and diverse molecules, and map the chemical space", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(21, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(22, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await session.step(23, "And the grid shows only the columns \"molregno, canonical_smiles\"", () => gridShowsColumns(page, "molregno, canonical_smiles"));
    await session.step(24, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(25, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(26, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(28, "When user picks \"Chem > Search > Similarity Search...\" from the top menu", () => pickFromTopMenu(page, "Chem > Search > Similarity Search..."), undefined, "Chem > Search > Similarity Search");
    await session.step(29, "Then Chem Similarity Search viewer should be visible", () => shouldBe(page, el("Chem Similarity Search viewer"), "visible"));
    await session.step(30, "When the video says \"The current molecule first, then the most similar, with scores\"", () => videoSays(page, "The current molecule first, then the most similar, with scores"));
    await session.step(31, "And user remembers the \"scores\" reading of Chem Similarity Search viewer", () => rememberReading(page, "scores", el("Chem Similarity Search viewer")));
    await session.step(33, "And user clicks on the \"cell 5 of canonical_smiles\" area of grid", () => clickArea(page, "cell 5 of canonical_smiles", el("grid")), undefined, "Make another molecule current: the search follows");
    await session.step(34, "Then row 5 should be current", () => currentRowIs(page, 5));
    await session.step(35, "And the \"scores\" reading of Chem Similarity Search viewer should not be as remembered", () => readingNotAsRemembered(page, "scores", el("Chem Similarity Search viewer")));
    await session.step(37, "When user clicks on card 2 of Chem Similarity Search viewer beside its drawing", () => clickNthCardBesideDrawing(page, 2, el("Chem Similarity Search viewer")), undefined, "Click a card to make its molecule current");
    await session.step(38, "Then the current row should be the row of the clicked card", () => currentRowIsClickedCard(page));
    await session.step(40, "When user picks \"Chem > Search > Diversity Search...\" from the top menu", () => pickFromTopMenu(page, "Chem > Search > Diversity Search..."), undefined, "Chem > Search > Diversity Search");
    await session.step(41, "Then Chem Diversity Search viewer should be visible", () => shouldBe(page, el("Chem Diversity Search viewer"), "visible"));
    await session.step(42, "When the video says \"Molecules picked to be as different from each other as possible\"", () => videoSays(page, "Molecules picked to be as different from each other as possible"));
    await session.step(44, "And user picks \"Chem > Analyze > Chemical Space...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Chemical Space..."), undefined, "Chem > Analyze > Chemical Space puts every molecule on a map");
    await session.step(45, "Then \"Chem Space\" dialog should be visible", () => shouldBe(page, el("\"Chem Space\" dialog"), "visible"));
    await session.step(46, "And Method input in \"Chem Space\" dialog should have value \"UMAP\"", () => shouldHaveValue(page, el("Method input in \"Chem Space\" dialog"), "UMAP"));
    await session.step(48, "When user clicks on OK button in \"Chem Space\" dialog", () => clickOn(page, el("OK button in \"Chem Space\" dialog")), undefined, "UMAP on fingerprint similarity, then OK");
    await session.step(49, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(50, "And a new column matching \"^Embed_X_\" should have been added", () => newColumnMatching(page, "^Embed_X_"));
    await session.step(51, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(52, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
    await session.step(53, "When the video says \"Each dot is a molecule; similar ones sit close together\"", () => videoSays(page, "Each dot is a molecule; similar ones sit close together"));
    await session.step(55, "And user hovers over the first free \"marker\" area of scatter plot viewer", () => hoverFirstFreeArea(page, "marker", el("scatter plot viewer")), undefined, "Point at a dot to see its structure");
    await session.step(56, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
    await session.step(57, "When user remembers the \"scores\" reading of Chem Similarity Search viewer", () => rememberReading(page, "scores", el("Chem Similarity Search viewer")));
    await session.step(59, "And user clicks on the \"marker of row 15\" area of scatter plot viewer", () => clickArea(page, "marker of row 15", el("scatter plot viewer")), undefined, "Click a dot: the similarity search follows it");
    await session.step(60, "Then row 15 should be current", () => currentRowIs(page, 15));
    await session.step(61, "And the \"scores\" reading of Chem Similarity Search viewer should not be as remembered", () => readingNotAsRemembered(page, "scores", el("Chem Similarity Search viewer")));
    await session.step(63, "When user picks \"NumAromaticCarbocycles\" in the \"color\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "NumAromaticCarbocycles", "color", el("scatter plot viewer")), undefined, "Color the map by any column");
    await session.step(64, "Then \"Color Column Name\" property of scatter plot viewer should be \"NumAromaticCarbocycles\"", () => propertyShouldBe(page, "Color Column Name", el("scatter plot viewer"), "NumAromaticCarbocycles"));
    await session.step(65, "And the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
  });
});
