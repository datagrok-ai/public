// import * as grok from 'datagrok-api/grok';
// import * as DG from 'datagrok-api/dg';

// let _progress: DG.TaskBarProgressIndicator | null = null;

// function geminiDownloadMonitor(m: CreateMonitor) {
//   m.addEventListener('downloadprogress', (e) => {
//     if (!_progress) {
//       grok.shell.info('Gemini model download started. This might take a few minutes.');
//       _progress = DG.TaskBarProgressIndicator.create('Downloading Gemini...');
//     }
//     _progress.update(e.loaded * 100, 'Downloading Gemini...');
//   });
// }

// export async function isHelpRequest(post: string): Promise<boolean> {
//   const schema = {
//     'type': 'boolean'
//   };

//   const systemPrompt = `You are a helpful LLM that inspects user input and returns true
//   when the input is likely a question concerning the functionality of the data analytics application.
//   Return false for anything else (random strings, identifiers, numbers, commands, requests)`;

//   const session = await LanguageModel.create({
//     monitor: geminiDownloadMonitor,
//     initialPrompts: [
//       {role: 'system', content: systemPrompt}
//     ]
//   });

//   const result = await session.prompt(post,
//     {
//       responseConstraint: schema,
//     }
//   );

//   return JSON.parse(result);
// }
