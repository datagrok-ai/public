import {_package} from './package';
import {KnimeHubClient} from './knime-hub-client';

let cachedClient: KnimeHubClient | null = null;

export function getKnimeClient(): KnimeHubClient {
  if (cachedClient)
    return cachedClient;

  const baseUrl = (_package.settings?.baseUrl ?? '') as string;
  if (!baseUrl)
    throw new Error('KNIME base URL is not configured in package settings');

  cachedClient = new KnimeHubClient(baseUrl);
  return cachedClient;
}
