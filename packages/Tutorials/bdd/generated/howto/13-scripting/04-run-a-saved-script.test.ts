/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/13-scripting/04-run-a-saved-script.feature
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
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {scriptPaneShows, scriptsSearchPutBack, variableReads} from '../../../bindings/howto-scripting.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, expand, hoverText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, scriptOnServer, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Scripting", () => {
  const session = feature(test, "features/howto/13-scripting/04-run-a-saved-script.feature", import.meta.url);
  test("Run a saved script without opening it", {tag: ["@guide", "@tutorial", "@viewport:1440x810"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And a script \"Body mass index\" is on the server:", () => scriptOnServer(page, "Body mass index", "//language: javascript\n//description: BMI from height and weight\n//input: double height = 175 {caption: Height, cm}\n//input: double weight = 70 {caption: Weight, kg}\n//output: double bmi\nbmi = weight / (height / 100) ** 2;"));
    await session.step(28, "And the Scripts gallery search is emptied at feature end", () => scriptsSearchPutBack(page));
    await session.step(29, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(30, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(33, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")), undefined, "Saved scripts are in Browse, under Platform > Functions > Scripts");
    await session.step(34, "And user expands \"Platform > Functions\" tree node inside browse tree", () => expand(page, el("\"Platform > Functions\" tree node inside browse tree")));
    await session.step(35, "And user clicks on \"Platform > Functions > Scripts\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Functions > Scripts\" tree node inside browse tree")));
    await session.step(36, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
    await session.step(37, "When user types \"Body mass\" into gallery search", () => typeInto(page, "Body mass", el("gallery search")));
    await session.step(39, "And user clicks on \"Body mass index\" link in gallery", () => clickOn(page, el("\"Body mass index\" link in gallery")), undefined, "Click a script: what it does, its inputs and outputs");
    await session.step(40, "Then the context panel should show \"Body mass index\"", () => contextPanelShows(page, "Body mass index"));
    await session.step(42, "When user clicks on \"Script\" pane header in context panel", () => clickOn(page, el("\"Script\" pane header in context panel")), undefined, "Script: its code, without opening an editor");
    await session.step(43, "Then the Script pane of the context panel should show \"bmi = weight / (height / 100) ** 2;\"", () => scriptPaneShows(page, "bmi = weight / (height / 100) ** 2;"));
    await session.step(45, "When user clicks on \"Run\" pane header in context panel", () => clickOn(page, el("\"Run\" pane header in context panel")), undefined, "Run: the same form, right in the context panel");
    await session.step(46, "And user enters \"160\" into \"Height, cm\" input in context panel", () => enterInto(page, "160", el("\"Height, cm\" input in context panel")));
    await session.step(47, "And user enters \"55\" into \"Weight, kg\" input in context panel", () => enterInto(page, "55", el("\"Weight, kg\" input in context panel")));
    await session.step(48, "And user clicks on context panel run button", () => clickOn(page, el("context panel run button")));
    await session.step(50, "Then the Variables panel should show \"bmi\" as \"21.48\"", () => variableReads(page, "bmi", "21.48"));
    await session.step(52, "When user hovers over the text \"bmi\" in \"Variables\" dock panel", () => hoverText(page, "bmi", el("\"Variables\" dock panel")), undefined, "The result goes to the Variables panel");
    await session.step(53, "Then no errors should have been logged", () => noErrors(page));
  });
});
