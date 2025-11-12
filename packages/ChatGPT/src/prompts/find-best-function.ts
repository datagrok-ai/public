import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';


/// Prompt API: https://developer.chrome.com/docs/ai/prompt-api
/// Enable all Gemini Nano-related flags here:
/// chrome://flags/#prompt-api-for-gemini-nano

let _progress: DG.TaskBarProgressIndicator | null = null;

export function geminiDownloadMonitor(m: CreateMonitor) {
  m.addEventListener('downloadprogress', (e) => {
    if (!_progress) {
      grok.shell.info('Gemini model download started. This might take a few minutes.');
      _progress = DG.TaskBarProgressIndicator.create('Downloading Gemini...');
    }
    _progress.update(e.loaded * 100, 'Downloading Gemini...');
  });
}

export async function findBestFunction(question: string): Promise<{function: string, parameterName: string} | null> {
  const schema = {
    "function": "string",
    "parameterName": "string"
  }
  
  const functions = [
    'Users joined after ${date}',
    'Bioactivity for bacterial targets for ${organism}',
  ]
  
  const systemPrompt = `You are a helpful LLM that matches a user question to one of the pre-defined 
  functions and extracts parameters, if possible. Parameters are specified as \${paramName}.
  You can only select at most one function. You should never return functions not in the list.
  You always return either null, or a valid JSON fitting the schema:' + JSONschema`;

  const session = await LanguageModel.create({
    monitor: geminiDownloadMonitor,
    initialPrompts: [
      { role: 'system', content: systemPrompt },
      { role: 'user', content: "Question: Shigella bioactivity\n\nFunctions:\n" + functions.join('\n') },
      { role: 'assistant', content: '{function: "Bioactivity for bacterial targets for ${organism}", parameter:"Shigella"}'},
    ]
  });

  const result = await session.prompt(
    `Question: ${question}\n\nFunctions:\n` + functions.join('\n'), 
    { responseConstraint: schema });
  
  return JSON.parse(result);
}
