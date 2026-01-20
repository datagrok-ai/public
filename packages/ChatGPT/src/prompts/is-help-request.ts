import {geminiDownloadMonitor} from '../llm-utils/query-matching';

export async function isHelpRequest(post: string): Promise<boolean> {
  const schema = {
    'type': 'boolean'
  };

  const systemPrompt = `You are a helpful LLM that inspects user input and returns true
  when the input is likely a question concerning the functionality of the data analytics application.
  Return false for anything else (random strings, identifiers, numbers, commands, requests)`;

  const session = await LanguageModel.create({
    monitor: geminiDownloadMonitor,
    initialPrompts: [
      {role: 'system', content: systemPrompt}
    ]
  });

  const result = await session.prompt(post,
    {
      responseConstraint: schema,
    }
  );

  return JSON.parse(result);
}
