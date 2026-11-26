/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/04-select-and-filter/05-hierarchical-filter.feature
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
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterIsExactlyCategory, filterPasses, filterPassesAll, noneOfFiltered, openFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickPanelMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Select and filter", () => {
  const session = feature(test, "features/howto/04-select-and-filter/05-hierarchical-filter.feature", import.meta.url);
  test("Filter with a hierarchy", {tag: ["@guide", "@tutorial", "@viewport:1280x720"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "And user opens the filter panel", () => openFilterPanel(page));
    await session.step(21, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(23, "When user picks \"Add Filter | Hierarchical\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Hierarchical"), undefined, "The panel's menu adds a hierarchical filter");
    await session.step(24, "Then \"SEX / DEMOG / USUBJID\" filter card should be visible", () => shouldBe(page, el("\"SEX / DEMOG / USUBJID\" filter card"), "visible"));
    await session.step(25, "When user hovers over \"SEX / DEMOG / USUBJID\" filter card", () => hoverOver(page, el("\"SEX / DEMOG / USUBJID\" filter card")));
    await session.step(27, "And user clicks on \"Select columns\" icon in \"SEX / DEMOG / USUBJID\" filter card", () => clickOn(page, el("\"Select columns\" icon in \"SEX / DEMOG / USUBJID\" filter card")), undefined, "The tree icon picks the columns, one per level");
    await session.step(28, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(29, "When user toggles the \"DEMOG\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "DEMOG", el("\"Select columns...\" dialog")));
    await session.step(30, "And user toggles the \"USUBJID\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "USUBJID", el("\"Select columns...\" dialog")));
    await session.step(31, "And user toggles the \"RACE\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "RACE", el("\"Select columns...\" dialog")));
    await session.step(32, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(33, "Then \"SEX / RACE\" filter card should be visible", () => shouldBe(page, el("\"SEX / RACE\" filter card"), "visible"));
    await session.step(34, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(36, "When user clicks on label of \"F\" hierarchy node", () => clickOn(page, el("label of \"F\" hierarchy node")), undefined, "Click a branch to keep only its rows");
    await session.step(37, "Then 3243 rows should pass the filter", () => filterPasses(page, 3243));
    await session.step(38, "And the filter should pass exactly the rows where \"SEX\" is \"F\"", () => filterIsExactlyCategory(page, "SEX", "F"));
    await session.step(40, "When user clicks on twistie of \"F\" hierarchy node", () => clickOn(page, el("twistie of \"F\" hierarchy node")), undefined, "The triangle opens the next level");
    await session.step(41, "Then \"Caucasian\" hierarchy node should be visible", () => shouldBe(page, el("\"Caucasian\" hierarchy node"), "visible"));
    await session.step(43, "When user clicks on checkbox of \"Caucasian\" hierarchy node", () => clickOn(page, el("checkbox of \"Caucasian\" hierarchy node")), undefined, "Untick a branch to leave it out");
    await session.step(44, "Then 420 rows should pass the filter", () => filterPasses(page, 420));
    await session.step(45, "And no rows where \"RACE\" is \"Caucasian\" should pass the filter", () => noneOfFiltered(page, "RACE", "Caucasian"));
    await session.step(46, "And no rows where \"SEX\" is \"M\" should pass the filter", () => noneOfFiltered(page, "SEX", "M"));
    await session.step(47, "And no errors should have been logged", () => noErrors(page));
  });
});
