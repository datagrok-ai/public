/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/markup/markup.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.markup]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, noErrors, propertyShouldBe, readingIs, readingReads, reportsNoError, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingExcludes, readingIncludes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Markup content and the interpretation mode the viewer resolved", () => {
  const session = feature(test, "features/viewers/markup/markup.feature", import.meta.url);
  test("Markup content and the interpretation mode the viewer resolved", {tag: ["@journey", "@viewers", "@realizes:viewers.markup"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(22, "And user adds a markup viewer", () => addViewer(page, "markup"));
    await session.step(23, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
    await session.step(24, "And the \"mode\" reading of markup viewer should be \"Markup\"", () => readingReads(page, "mode", el("markup viewer"), "Markup"));
    await session.step(25, "And \"mode\" property of markup viewer should be \"Auto\"", () => propertyShouldBe(page, "mode", el("markup viewer"), "Auto"));
    await session.step(26, "And the \"markup enabled\" reading of markup viewer should be \"true\"", () => readingReads(page, "markup enabled", el("markup viewer"), "true"));
    await session.step(27, "And markup viewer should report no error", () => reportsNoError(page, el("markup viewer")));
    await run.scenario("The viewer opens on the Markdown sample rendered, not on its source", async () => {
      await session.step(30, "Then the \"heading 1\" reading of markup viewer should be \"What’s Markdown?\"", () => readingReads(page, "heading 1", el("markup viewer"), "What’s Markdown?"));
      await session.step(31, "And the \"list items\" reading of markup viewer should be 4", () => readingIs(page, "list items", el("markup viewer"), 4));
      await session.step(32, "And the \"links\" reading of markup viewer should be 4", () => readingIs(page, "links", el("markup viewer"), 4));
      await session.step(33, "And the \"preformatted\" reading of markup viewer should be 0", () => readingIs(page, "preformatted", el("markup viewer"), 0));
      await session.step(34, "And the \"text\" reading of markup viewer should include the text \"Markdown is a lightweight markup language\"", () => readingIncludes(page, "text", el("markup viewer"), "Markdown is a lightweight markup language"));
      await session.step(35, "And the \"text\" reading of markup viewer should not include the text \"#\"", () => readingExcludes(page, "text", el("markup viewer"), "#"));
      await session.step(36, "And the \"text\" reading of markup viewer should not include the text \"https://\"", () => readingExcludes(page, "text", el("markup viewer"), "https://"));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Auto resolves the mode from the content, and the property never says which", async () => {
      await session.step(40, "When user sets \"content\" property of markup viewer to \"# Heading probe\"", () => setProperty(page, "content", el("markup viewer"), "# Heading probe"));
      await session.step(41, "Then \"mode\" property of markup viewer should be \"Auto\"", () => propertyShouldBe(page, "mode", el("markup viewer"), "Auto"));
      await session.step(42, "And the \"mode\" reading of markup viewer should be \"Markup\"", () => readingReads(page, "mode", el("markup viewer"), "Markup"));
      await session.step(43, "And the \"heading 1\" reading of markup viewer should be \"Heading probe\"", () => readingReads(page, "heading 1", el("markup viewer"), "Heading probe"));
      await session.step(44, "When user sets \"content\" property of markup viewer to \"<b>bold probe</b> plain probe\"", () => setProperty(page, "content", el("markup viewer"), "<b>bold probe</b> plain probe"));
      await session.step(45, "Then \"mode\" property of markup viewer should be \"Auto\"", () => propertyShouldBe(page, "mode", el("markup viewer"), "Auto"));
      await session.step(46, "And the \"mode\" reading of markup viewer should be \"Html\"", () => readingReads(page, "mode", el("markup viewer"), "Html"));
      await session.step(47, "And the \"bold\" reading of markup viewer should be 1", () => readingIs(page, "bold", el("markup viewer"), 1));
      await session.step(48, "And the \"heading 1\" reading of markup viewer should be \"\"", () => readingReads(page, "heading 1", el("markup viewer"), ""));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Mode None shows the viewer's own markup verbatim, escaped, in a pre", async () => {
      await session.step(52, "When user sets \"content\" property of markup viewer to \"<b>bold probe</b> plain probe\"", () => setProperty(page, "content", el("markup viewer"), "<b>bold probe</b> plain probe"));
      await session.step(53, "And user sets \"mode\" property of markup viewer to \"None\"", () => setProperty(page, "mode", el("markup viewer"), "None"));
      await session.step(54, "Then the \"mode\" reading of markup viewer should be \"None\"", () => readingReads(page, "mode", el("markup viewer"), "None"));
      await session.step(55, "And the \"preformatted\" reading of markup viewer should be 1", () => readingIs(page, "preformatted", el("markup viewer"), 1));
      await session.step(56, "And the \"bold\" reading of markup viewer should be 0", () => readingIs(page, "bold", el("markup viewer"), 0));
      await session.step(57, "And the \"text\" reading of markup viewer should include the text \"<b>bold probe</b> plain probe\"", () => readingIncludes(page, "text", el("markup viewer"), "<b>bold probe</b> plain probe"));
      await session.step(58, "When user sets \"mode\" property of markup viewer to \"Auto\"", () => setProperty(page, "mode", el("markup viewer"), "Auto"));
      await session.step(59, "Then the \"mode\" reading of markup viewer should be \"Html\"", () => readingReads(page, "mode", el("markup viewer"), "Html"));
      await session.step(60, "And the \"preformatted\" reading of markup viewer should be 0", () => readingIs(page, "preformatted", el("markup viewer"), 0));
      await session.step(61, "And the \"bold\" reading of markup viewer should be 1", () => readingIs(page, "bold", el("markup viewer"), 1));
      await session.step(62, "And the \"text\" reading of markup viewer should not include the text \"<b>\"", () => readingExcludes(page, "text", el("markup viewer"), "<b>"));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Mode Html leaves a hash a hash; Mode Markup makes it a heading whatever the content looks like", async () => {
      await session.step(66, "When user sets \"content\" property of markup viewer to \"# Heading probe\"", () => setProperty(page, "content", el("markup viewer"), "# Heading probe"));
      await session.step(67, "And user sets \"mode\" property of markup viewer to \"Html\"", () => setProperty(page, "mode", el("markup viewer"), "Html"));
      await session.step(68, "Then the \"mode\" reading of markup viewer should be \"Html\"", () => readingReads(page, "mode", el("markup viewer"), "Html"));
      await session.step(69, "And the \"heading 1\" reading of markup viewer should be \"\"", () => readingReads(page, "heading 1", el("markup viewer"), ""));
      await session.step(70, "And the \"text\" reading of markup viewer should include the text \"# Heading probe\"", () => readingIncludes(page, "text", el("markup viewer"), "# Heading probe"));
      await session.step(71, "When user sets \"mode\" property of markup viewer to \"Markup\"", () => setProperty(page, "mode", el("markup viewer"), "Markup"));
      await session.step(72, "Then the \"mode\" reading of markup viewer should be \"Markup\"", () => readingReads(page, "mode", el("markup viewer"), "Markup"));
      await session.step(73, "And the \"heading 1\" reading of markup viewer should be \"Heading probe\"", () => readingReads(page, "heading 1", el("markup viewer"), "Heading probe"));
      await session.step(74, "And the \"text\" reading of markup viewer should not include the text \"#\"", () => readingExcludes(page, "text", el("markup viewer"), "#"));
      await session.step(75, "When user sets \"mode\" property of markup viewer to \"Auto\"", () => setProperty(page, "mode", el("markup viewer"), "Auto"));
      await session.step(76, "Then the \"mode\" reading of markup viewer should be \"Markup\"", () => readingReads(page, "mode", el("markup viewer"), "Markup"));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The markup pass turns the emphasis marks into a strong element and eats them", async () => {
      await session.step(80, "When user sets \"mode\" property of markup viewer to \"Markup\"", () => setProperty(page, "mode", el("markup viewer"), "Markup"));
      await session.step(81, "And user sets \"content\" property of markup viewer to \"**md bold** plain tail\"", () => setProperty(page, "content", el("markup viewer"), "**md bold** plain tail"));
      await session.step(82, "Then the \"bold\" reading of markup viewer should be 1", () => readingIs(page, "bold", el("markup viewer"), 1));
      await session.step(83, "And the \"text\" reading of markup viewer should be \"md bold plain tail\"", () => readingReads(page, "text", el("markup viewer"), "md bold plain tail"));
      await session.step(84, "And the \"text\" reading of markup viewer should not include the text \"**\"", () => readingExcludes(page, "text", el("markup viewer"), "**"));
      await session.step(85, "When user sets \"mode\" property of markup viewer to \"Auto\"", () => setProperty(page, "mode", el("markup viewer"), "Auto"));
      await session.step(86, "Then the \"bold\" reading of markup viewer should be 1", () => readingIs(page, "bold", el("markup viewer"), 1));
      await session.step(87, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Markdown list is a list of items, not four lines of source", async () => {
      await session.step(90, "When user sets \"content\" property of markup viewer to \"* one\\n* two\\n* three\"", () => setProperty(page, "content", el("markup viewer"), "* one\\n* two\\n* three"));
      await session.step(91, "Then the \"mode\" reading of markup viewer should be \"Markup\"", () => readingReads(page, "mode", el("markup viewer"), "Markup"));
      await session.step(92, "And the \"list items\" reading of markup viewer should be 3", () => readingIs(page, "list items", el("markup viewer"), 3));
      await session.step(93, "And the \"preformatted\" reading of markup viewer should be 0", () => readingIs(page, "preformatted", el("markup viewer"), 0));
      await session.step(94, "And the \"text\" reading of markup viewer should not include the text \"*\"", () => readingExcludes(page, "text", el("markup viewer"), "*"));
      await session.step(95, "When user sets \"mode\" property of markup viewer to \"None\"", () => setProperty(page, "mode", el("markup viewer"), "None"));
      await session.step(96, "Then the \"list items\" reading of markup viewer should be 0", () => readingIs(page, "list items", el("markup viewer"), 0));
      await session.step(97, "And the \"preformatted\" reading of markup viewer should be 1", () => readingIs(page, "preformatted", el("markup viewer"), 1));
      await session.step(98, "And the \"text\" reading of markup viewer should include the text \"* one\"", () => readingIncludes(page, "text", el("markup viewer"), "* one"));
      await session.step(99, "When user sets \"mode\" property of markup viewer to \"Auto\"", () => setProperty(page, "mode", el("markup viewer"), "Auto"));
      await session.step(100, "Then the \"list items\" reading of markup viewer should be 3", () => readingIs(page, "list items", el("markup viewer"), 3));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Markup Enabled decides whether the table expressions are evaluated at all", async () => {
      await session.step(104, "When user sets \"content\" property of markup viewer to \"Rows: #{t.rowCount}\"", () => setProperty(page, "content", el("markup viewer"), "Rows: #{t.rowCount}"));
      await session.step(105, "Then the \"text\" reading of markup viewer should be \"Rows: 1000\"", () => readingReads(page, "text", el("markup viewer"), "Rows: 1000"));
      await session.step(106, "When user sets \"markupEnabled\" property of markup viewer to \"false\"", () => setProperty(page, "markupEnabled", el("markup viewer"), "false"));
      await session.step(107, "Then the \"markup enabled\" reading of markup viewer should be \"false\"", () => readingReads(page, "markup enabled", el("markup viewer"), "false"));
      await session.step(108, "And the \"text\" reading of markup viewer should be \"Rows: #{t.rowCount}\"", () => readingReads(page, "text", el("markup viewer"), "Rows: #{t.rowCount}"));
      await session.step(109, "When user sets \"markupEnabled\" property of markup viewer to \"true\"", () => setProperty(page, "markupEnabled", el("markup viewer"), "true"));
      await session.step(110, "Then the \"markup enabled\" reading of markup viewer should be \"true\"", () => readingReads(page, "markup enabled", el("markup viewer"), "true"));
      await session.step(111, "And the \"text\" reading of markup viewer should be \"Rows: 1000\"", () => readingReads(page, "text", el("markup viewer"), "Rows: 1000"));
      await session.step(112, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
