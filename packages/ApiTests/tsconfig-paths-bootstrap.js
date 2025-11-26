import tsConfigPaths from 'tsconfig-paths';
import { readFileSync } from 'fs';
import path from 'path';

const tsConfig = JSON.parse(readFileSync('./tsconfig.node.json', 'utf8'));
const baseUrl = path.resolve('./dist-node');
const paths = tsConfig.compilerOptions.paths;

tsConfigPaths.register({ baseUrl, paths });
