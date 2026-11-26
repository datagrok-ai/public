/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/13-scripting/05-learn-from-the-samples.feature
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
import {clickOn, expand, hoverOver, hoverText, pressKey, replaceLine, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {colorCodedAs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, contextPanelClosed, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Scripting", () => {
  const session = feature(test, "features/howto/13-scripting/05-learn-from-the-samples.feature", import.meta.url);
  test("Learn from the sample scripts", {tag: ["@guide", "@tutorial", "@viewport:1440x810"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")), undefined, "JS API Samples: short working scripts, by topic");
    await session.step(22, "And user expands \"Platform > Functions\" tree node inside browse tree", () => expand(page, el("\"Platform > Functions\" tree node inside browse tree")));
    await session.step(23, "And user expands \"Platform > Functions > JS API Samples\" tree node inside browse tree", () => expand(page, el("\"Platform > Functions > JS API Samples\" tree node inside browse tree")));
    await session.step(24, "And user expands \"Platform > Functions > JS API Samples > grid\" tree node inside browse tree", () => expand(page, el("\"Platform > Functions > JS API Samples > grid\" tree node inside browse tree")));
    await session.step(25, "And user expands \"Platform > Functions > JS API Samples > grid > color-coding\" tree node inside browse tree", () => expand(page, el("\"Platform > Functions > JS API Samples > grid > color-coding\" tree node inside browse tree")));
    await session.step(27, "And user clicks on \"Platform > Functions > JS API Samples > grid > color-coding > color-coding\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Functions > JS API Samples > grid > color-coding > color-coding\" tree node inside browse tree")), undefined, "A sample opens as a script you can change and run");
    await session.step(28, "Then code editor should contain the text \"t.col('height').meta.colors.setConditional\"", () => shouldContainText(page, el("code editor"), "t.col('height').meta.colors.setConditional"));
    await session.step(30, "When user hovers over the text \"// Conditional color-coding for numerical columns\" in code editor", () => hoverText(page, "// Conditional color-coding for numerical columns", el("code editor")), undefined, "It colors the columns of a demo table, one line per column");
    await session.step(32, "And user replaces the line starting with \"t.col('height')\" in code editor with \"t.col('height').meta.colors.setConditional({'20-160': '#00FF00', '160-190': '#220505'});\"", () => replaceLine(page, "t.col('height')", el("code editor"), "t.col('height').meta.colors.setConditional({'20-160': '#00FF00', '160-190': '#220505'});"), undefined, "Change it: green up to 160 cm, dark from 160 to 190");
    await session.step(34, "And user hovers over \"Run script (F5)\" icon", () => hoverOver(page, el("\"Run script (F5)\" icon")), undefined, "Run it: this play icon, or F5");
    await session.step(36, "And user hovers over the text \"addTableView\" in code editor", () => hoverText(page, "addTableView", el("code editor")), undefined, "Its last line opens the table in a view");
    await session.step(37, "And user presses F5", () => pressKey(page, "F5"));
    await session.step(38, "Then the \"demog 10000\" view should be current", () => viewIsCurrent(page, "demog 10000"));
    await session.step(39, "And \"height\" column should be color-coded conditionally", () => colorCodedAs(page, "height", "conditionally"));
    await session.step(40, "And \"age\" column should be color-coded linearly", () => colorCodedAs(page, "age", "linearly"));
    await session.step(42, "When user hovers over the \"header height\" area of grid", () => hoverArea(page, "header height", el("grid")), undefined, "The script opened a demo table with its colors");
    await session.step(43, "Then no errors should have been logged", () => noErrors(page));
  });
});
