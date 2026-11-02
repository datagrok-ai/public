/** The AI features, as the prefixes of their session IDs: the usage log reads a turn's feature from its session. */
export enum AiFeature {
  Chat = 'chat',
  Engine = 'engine',
  Execute = 'execute',
  Help = 'help',
  Search = 'search',
  QueryMatch = 'query-match',
  DbIndex = 'db-index',
  Benchmark = 'benchmark',
}

/** Chats restored from before these prefixes start with `claude-`, so anything unmatched is a chat. */
export function featureOf(sessionId: string): AiFeature {
  return Object.values(AiFeature).find((f) => sessionId.startsWith(`${f}-`)) ?? AiFeature.Chat;
}
