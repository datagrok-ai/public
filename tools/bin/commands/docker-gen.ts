import {generateCeleryArtifacts} from '../utils/python-celery-gen';
import * as color from '../utils/color-utils';

export async function dockerGen(args: any) {
  color.setVerbose(args.verbose || args.v || false);
  const packageDir = process.cwd();

  const generated = generateCeleryArtifacts(packageDir);
  if (generated)
    color.success('Generated Celery Docker artifacts');
  else
    console.log('No python/ directory with annotated functions found');

  return true;
}
