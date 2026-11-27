/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/08-end-to-end/03-compare-groups.feature
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
import {ownLayoutsGone} from '../../../bindings/howto-dashboards.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, hoverOver, pressKey, selectIn, shouldBe, shouldContainText, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterIsExactlyCategory} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, layoutsDeleted, openDatasetRowsAs, simpleModeOff, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hasArea, hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · End to end", () => {
  const session = feature(test, "features/howto/08-end-to-end/03-compare-groups.feature", import.meta.url);
  test("Compare groups", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And simple mode is off", () => simpleModeOff(page));
    await session.step(22, "And the layouts named \"Demographics, Weight by sex and race\" are deleted when the feature ends", () => layoutsDeleted(page, "Demographics, Weight by sex and race"));
    await session.step(23, "And the user's own layouts named \"Demographics, Weight by sex and race\" are removed now and at feature end", () => ownLayoutsGone(page, "Demographics, Weight by sex and race"));
    await session.step(24, "And user opens demog dataset keeping the first 5850 rows as \"Demographics\"", () => openDatasetRowsAs(page, ds("demog"), 5850, "Demographics"));
    await session.step(25, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(26, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(27, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(28, "When the video says \"The question: do men and women differ in weight, in every race\"", () => videoSays(page, "The question: do men and women differ in weight, in every race"));
    await session.step(30, "And user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")), undefined, "First choose who to compare: the RA patients");
    await session.step(31, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
    await session.step(32, "When user clicks on the \"category RA of DIS_POP\" area of filter panel", () => clickArea(page, "category RA of DIS_POP", el("filter panel")));
    await session.step(33, "Then the filter should pass exactly the rows where \"DIS_POP\" is \"RA\"", () => filterIsExactlyCategory(page, "DIS_POP", "RA"));
    await session.step(35, "When user clicks on box plot icon on toolbox", () => clickOn(page, el("box plot icon on toolbox")), undefined, "A box plot compares weight between the sexes");
    await session.step(36, "Then box plot viewer should be visible", () => shouldBe(page, el("box plot viewer"), "visible"));
    await session.step(37, "When user picks \"WEIGHT\" in the \"Value\" column selector of box plot viewer", () => pickInColumnSelector(page, "WEIGHT", "Value", el("box plot viewer")));
    await session.step(38, "And user picks \"SEX\" in the \"Category 1\" column selector of box plot viewer", () => pickInColumnSelector(page, "SEX", "Category 1", el("box plot viewer")));
    await session.step(40, "And user hovers over the \"p value\" area of box plot viewer", () => hoverArea(page, "p value", el("box plot viewer")), undefined, "Men weigh more; the tiny p-value says it is unlikely to be chance");
    await session.step(41, "Then tooltip should contain text \"Welch's t-test\"", () => shouldContainText(page, el("tooltip"), "Welch's t-test"));
    await session.step(43, "When user picks \"RACE\" in the \"Category 1\" column selector of box plot viewer", () => pickInColumnSelector(page, "RACE", "Category 1", el("box plot viewer")), undefined, "Now check every race: make race the category");
    await session.step(45, "And user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")), undefined, "and split each race by sex, in the settings");
    await session.step(46, "And user clicks on settings icon of box plot viewer", () => clickOn(page, el("settings icon of box plot viewer")));
    await session.step(47, "Then the context panel should show \"Box plot\"", () => contextPanelShows(page, "Box plot"));
    await session.step(48, "When user selects \"SEX\" in \"Category2\" property in context panel", () => selectIn(page, "SEX", el("\"Category2\" property in context panel")));
    await session.step(49, "Then box plot viewer should have a \"category Asian | F\" area", () => hasArea(page, el("box plot viewer"), "category Asian | F"));
    await session.step(50, "When the video says \"Within every race, the men's box sits higher\"", () => videoSays(page, "Within every race, the men's box sits higher"));
    await session.step(51, "And user hovers over the \"category Black | M\" area of box plot viewer", () => hoverArea(page, "category Black | M", el("box plot viewer")));
    await session.step(53, "And user expands Layouts accordion header in toolbox", () => expand(page, el("Layouts accordion header in toolbox")), undefined, "Keep this arrangement: Toolbox > Layouts > Save");
    await session.step(54, "And user clicks on Save button in layouts pane", () => clickOn(page, el("Save button in layouts pane")));
    await session.step(55, "Then \"Demographics\" layout card should be visible", () => shouldBe(page, el("\"Demographics\" layout card"), "visible"));
    await session.step(56, "And the context panel should show \"Demographics\"", () => contextPanelShows(page, "Demographics"));
    await session.step(58, "When user types \"Weight by sex and race\" into Name details field in context panel", () => typeInto(page, "Weight by sex and race", el("Name details field in context panel")), undefined, "Give it a name you will find again");
    await session.step(59, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(60, "Then \"Weight by sex and race\" layout card should be visible", () => shouldBe(page, el("\"Weight by sex and race\" layout card"), "visible"));
    await session.step(61, "When the video says \"Its card applies it to any table with these columns\"", () => videoSays(page, "Its card applies it to any table with these columns"));
    await session.step(62, "And user hovers over \"Weight by sex and race\" layout card", () => hoverOver(page, el("\"Weight by sex and race\" layout card")));
    await session.step(63, "Then no errors should have been logged", () => noErrors(page));
  });
});
