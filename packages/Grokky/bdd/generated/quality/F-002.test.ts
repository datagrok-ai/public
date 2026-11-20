/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/quality/F-002.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {aiPanelNewChat, askAI, liveModelGate, runtimeHealthy, screenshot, trafficRecorded, treeMapColoredBy} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A tree map asked to color by a text column shows its categories instead of grey", () => {
  const session = feature(test, "features/quality/F-002.feature", import.meta.url);
  test("The tree map asked to color by disease puts disease last in its split and sets no color column", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(13, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(14, "And user is logged in", () => loggedIn(page));
    await session.step(15, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(16, "And user opens a table \"clinic\" with:", () => openTableOf(page, "clinic", [["id","sex","site","disease","age"],["1","Female","Los Angeles","Anxiety","34"],["2","Female","Los Angeles","Anxiety","51"],["3","Female","Los Angeles","Anxiety","47"],["4","Female","Los Angeles","Anxiety","62"],["5","Female","Los Angeles","Anxiety","29"],["6","Female","Buffalo","Anxiety","58"],["7","Female","Buffalo","Anxiety","41"],["8","Female","Buffalo","Anxiety","66"],["9","Male","Los Angeles","Anxiety","38"],["10","Male","Los Angeles","Anxiety","55"],["11","Male","New York","Glaucoma","71"],["12","Male","New York","Glaucoma","49"],["13","Male","New York","Glaucoma","63"],["14","Male","New York","Glaucoma","57"],["15","Male","New York","Hepatitis A","44"],["16","Male","New York","Hepatitis A","36"],["17","Male","New York","Hepatitis A","52"],["18","Male","Orlando","Hepatitis A","68"],["19","Male","Orlando","Hepatitis A","45"],["20","Male","Orlando","Hepatitis A","59"],["21","Male","Orlando","Hepatitis A","33"],["22","Male","Buffalo","Anxiety","61"],["23","Male","Buffalo","Anxiety","50"]]), [["id","sex","site","disease","age"],["1","Female","Los Angeles","Anxiety","34"],["2","Female","Los Angeles","Anxiety","51"],["3","Female","Los Angeles","Anxiety","47"],["4","Female","Los Angeles","Anxiety","62"],["5","Female","Los Angeles","Anxiety","29"],["6","Female","Buffalo","Anxiety","58"],["7","Female","Buffalo","Anxiety","41"],["8","Female","Buffalo","Anxiety","66"],["9","Male","Los Angeles","Anxiety","38"],["10","Male","Los Angeles","Anxiety","55"],["11","Male","New York","Glaucoma","71"],["12","Male","New York","Glaucoma","49"],["13","Male","New York","Glaucoma","63"],["14","Male","New York","Glaucoma","57"],["15","Male","New York","Hepatitis A","44"],["16","Male","New York","Hepatitis A","36"],["17","Male","New York","Hepatitis A","52"],["18","Male","Orlando","Hepatitis A","68"],["19","Male","Orlando","Hepatitis A","45"],["20","Male","Orlando","Hepatitis A","59"],["21","Male","Orlando","Hepatitis A","33"],["22","Male","Buffalo","Anxiety","61"],["23","Male","Buffalo","Anxiety","50"]]);
    await session.step(41, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(44, "When user asks the AI panel \"say I just got this clinical table from a colleague. What is in it?\"", () => askAI(page, "say I just got this clinical table from a colleague. What is in it?"));
    await session.step(45, "And user asks the AI panel \"Thanks. What's worth a closer look? Anything odd or suspicious in the data?\"", () => askAI(page, "Thanks. What's worth a closer look? Anything odd or suspicious in the data?"));
    await session.step(46, "And user asks the AI panel \"ok, show me the sex / site / disease overlap in a tree map colored by disease, and tell me what it shows\"", () => askAI(page, "ok, show me the sex / site / disease overlap in a tree map colored by disease, and tell me what it shows"));
    await session.step(47, "Then the open tableview should have 1 Tree map viewer", () => viewerCount(page, 1, "Tree map"));
    await session.step(48, "And the tree map should color its blocks by \"disease\"", () => treeMapColoredBy(page, "disease"));
    await session.step(49, "And a screenshot \"F-002-tree-map\" of the AI panel is saved", () => screenshot(page, "F-002-tree-map"));
  });
});
