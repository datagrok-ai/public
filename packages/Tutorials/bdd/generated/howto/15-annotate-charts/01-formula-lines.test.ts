/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/15-annotate-charts/01-formula-lines.feature
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
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {viewerShare} from '../../../bindings/howto-viewers.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {drawsFormulaLines, formulaLineRanges} from '@datagrok-libraries/bdd/bindings/tiers/viewers/formula-lines';
import {addViewerWith, areaColor, noErrors, pickFromAreaContextMenu, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Annotations", () => {
  const session = feature(test, "features/howto/15-annotate-charts/01-formula-lines.feature", import.meta.url);
  test("Draw reference lines with formula lines", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","HEIGHT"],["Y","WEIGHT"]]), [["X","HEIGHT"],["Y","WEIGHT"]]);
    await session.step(23, "And scatter plot viewer takes 65% of the space it shares with the grid", () => viewerShare(page, el("scatter plot viewer"), 65));
    await session.step(24, "When the video says \"Formula lines mark a value or a curve on a chart\"", () => videoSays(page, "Formula lines mark a value or a curve on a chart"));
    await session.step(26, "And user picks \"Annotations > Add Line\" from the context menu of the \"y axis\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Annotations > Add Line", "y axis", el("scatter plot viewer")), undefined, "Right-click an axis: Annotations > Add Line");
    await session.step(27, "Then \"Formula Lines\" dialog should be visible", () => shouldBe(page, el("\"Formula Lines\" dialog"), "visible"));
    await session.step(29, "When user enters \"100\" into \"Value\" text input in \"Formula Lines\" dialog", () => enterInto(page, "100", el("\"Value\" text input in \"Formula Lines\" dialog")), undefined, "Type the value the line marks");
    await session.step(31, "And user clicks on \"Add new\" button in \"Formula Lines\" dialog", () => clickOn(page, el("\"Add new\" button in \"Formula Lines\" dialog")), undefined, "Add new > Line takes a formula of X and Y");
    await session.step(32, "And user picks \"Line\" from the open menu", () => pickFromOpenMenu(page, "Line"));
    await session.step(34, "And user enters \"${WEIGHT} = 25 * (${HEIGHT} / 100) ^ 2\" into text area in \"Line\" section", () => enterInto(page, "${WEIGHT} = 25 * (${HEIGHT} / 100) ^ 2", el("text area in \"Line\" section")), undefined, "Here: the weight at a BMI of 25");
    await session.step(36, "And user enters \"BMI 25\" into \"Title\" text input in \"Formula Lines\" dialog", () => enterInto(page, "BMI 25", el("\"Title\" text input in \"Formula Lines\" dialog")), undefined, "The title is written on the line");
    await session.step(38, "And user enters \"#d62728\" into \"Color\" input in \"Formula Lines\" dialog", () => enterInto(page, "#d62728", el("\"Color\" input in \"Formula Lines\" dialog")), undefined, "Format: the line's color, style and width");
    await session.step(39, "And user clicks on OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(41, "Then scatter plot viewer should draw 2 formula lines", () => drawsFormulaLines(page, el("scatter plot viewer"), 2), undefined, "Points above the curve have a BMI over 25");
    await session.step(42, "And the \"view\" area of scatter plot viewer should contain the color \"#d62728\"", () => areaColor(page, "view", el("scatter plot viewer"), "#d62728"));
    await session.step(44, "When user picks \"Tools > Formula Lines...\" from the context menu of the \"view\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Tools > Formula Lines...", "view", el("scatter plot viewer")), undefined, "Tools > Formula Lines... lists and edits them all");
    await session.step(45, "Then \"Formula Lines\" dialog should be visible", () => shouldBe(page, el("\"Formula Lines\" dialog"), "visible"));
    await session.step(46, "When user clicks on CANCEL button in \"Formula Lines\" dialog", () => clickOn(page, el("CANCEL button in \"Formula Lines\" dialog")));
    await session.step(47, "Then the formula lines of scatter plot viewer should hold \"${WEIGHT} = 100\" over the ranges \"..\"", () => formulaLineRanges(page, el("scatter plot viewer"), "${WEIGHT} = 100", ".."));
    await session.step(48, "And the formula lines of scatter plot viewer should hold \"${WEIGHT} = 25 * (${HEIGHT} / 100) ^ 2\" over the ranges \"..\"", () => formulaLineRanges(page, el("scatter plot viewer"), "${WEIGHT} = 25 * (${HEIGHT} / 100) ^ 2", ".."));
    await session.step(49, "And no errors should have been logged", () => noErrors(page));
  });
});
