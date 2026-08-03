import {NodeApiClient} from './node-dapi';
import {getDevKey} from './test-utils';

export async function createClient(hostArg?: string): Promise<NodeApiClient> {
  const {url, key} = getDevKey(hostArg ?? '');
  return NodeApiClient.login(url, key);
}
