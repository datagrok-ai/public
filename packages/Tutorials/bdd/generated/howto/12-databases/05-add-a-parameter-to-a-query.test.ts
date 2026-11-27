/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/12-databases/05-add-a-parameter-to-a-query.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
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
import {runWithPlayIcon} from '../../../bindings/howto-databases.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {appendToEditor, clickOn, enterInto, expand, hoverOver, hoverText, replaceCode, shouldBe, shouldHaveValue, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelClosed, currentViewType, dialogCloses, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Databases", () => {
  const session = feature(test, "features/howto/12-databases/05-add-a-parameter-to-a-query.feature", import.meta.url);
  test("Add a parameter to a query", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(16, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "When the video says \"A parameter makes one query answer many questions\"", () => videoSays(page, "A parameter makes one query answer many questions"));
    await session.step(19, "And user expands Databases tree node inside browse tree", () => expand(page, el("Databases tree node inside browse tree")));
    await session.step(20, "And user expands Databases---PostgresDart tree node inside browse tree", () => expand(page, el("Databases---PostgresDart tree node inside browse tree")));
    await session.step(22, "And user picks \"New Query...\" from the context menu of Databases---PostgresDart---Biologics tree node inside browse tree", () => pickFromContextMenu(page, "New Query...", el("Databases---PostgresDart---Biologics tree node inside browse tree")), undefined, "A new query on the connection");
    await session.step(23, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
    await session.step(24, "When user enters \"Assay types in a category\" into Name input", () => enterInto(page, "Assay types in a category", el("Name input")));
    await session.step(25, "And user replaces the code of code editor with '--input: string category = \"Developability\"'", () => replaceCode(page, el("code editor"), "--input: string category = \"Developability\""));
    await session.step(27, "And user hovers over the text \"--input\" in code editor", () => hoverText(page, "--input", el("code editor")), undefined, "An --input line declares the parameter and its default");
    await session.step(28, "And user appends \"select name, description from biologics.assay_types where category = @category\" to code editor", () => appendToEditor(page, "select name, description from biologics.assay_types where category = @category", el("code editor")));
    await session.step(30, "And user hovers over the text \"@category\" in code editor", () => hoverText(page, "@category", el("code editor")), undefined, "The SQL uses it as @category");
    await session.step(32, "And user runs the query with the play icon", () => runWithPlayIcon(page), undefined, "Run it: it asks for the value");
    await session.step(33, "Then \"Assay types in a category\" dialog should be visible", () => shouldBe(page, el("\"Assay types in a category\" dialog"), "visible"));
    await session.step(34, "And Category input in \"Assay types in a category\" dialog should have value \"Developability\"", () => shouldHaveValue(page, el("Category input in \"Assay types in a category\" dialog"), "Developability"));
    await session.step(36, "When user hovers over Category input in \"Assay types in a category\" dialog", () => hoverOver(page, el("Category input in \"Assay types in a category\" dialog")), undefined, "It starts with the default, Developability");
    await session.step(38, "And user enters \"Legacy\" into Category input in \"Assay types in a category\" dialog", () => enterInto(page, "Legacy", el("Category input in \"Assay types in a category\" dialog")), undefined, "Type another category");
    await session.step(39, "And user clicks on OK button in \"Assay types in a category\" dialog", () => clickOn(page, el("OK button in \"Assay types in a category\" dialog")));
    await session.step(40, "Then the \"Assay types in a category\" dialog should close", () => dialogCloses(page, "Assay types in a category"));
    await session.step(41, "And grid should be visible", () => shouldBe(page, el("grid"), "visible"));
    await session.step(42, "And the \"rows\" reading of grid should be 9", () => readingIs(page, "rows", el("grid"), 9));
    await session.step(43, "When the video says \"Only the Legacy assay types come back\" at the bottom", () => videoSaysAt(page, "Only the Legacy assay types come back", "bottom"));
    await session.step(44, "And user hovers over grid", () => hoverOver(page, el("grid")));
    await session.step(45, "Then no errors should have been logged", () => noErrors(page));
  });
});
