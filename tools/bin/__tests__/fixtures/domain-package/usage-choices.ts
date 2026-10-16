// Positive tsc fixture: a choices alias holds the manifest's values verbatim, quotes and backslashes included.
import {SampleResolution} from './src/generated/db';

export const values: SampleResolution[] = ['Fixed', `Won't fix`, 'a\\b'];
