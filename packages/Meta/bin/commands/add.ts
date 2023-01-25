import fs from 'fs';
import path from 'path';

function addPackageVersion(name: string, description: string, packageVersion: string, dependency: string,
                           dependencyVersion: string, repository: Map, category: string, rebuild: boolean) {

  const jsonTemplateLoc = path.join(path.dirname(path.dirname(__dirname)), 'package.json');
  let jsonContent = JSON.parse(fs.readFileSync(jsonTemplateLoc, 'utf8'));

  var data = {
    [name]: {
      dependencies: {
        [packageVersion]: {}
      },
      description: description
    }
  }

  if (repository) {
    data[name]['repository'] = repository;
  } else {
    var err = 'The repository should be specified';
    if (jsonContent.data.hasOwnProperty(name)) {
      var presult = jsonContent['data'][name]
      if (!presult.hasOwnProperty('repository')) {
        throw err;
      }
    } else {
      throw err;
    }
  }

  if (category) {
    data[name]['category'] = category;
  }

  if (dependency && dependencyVersion) {
    data[name]['dependencies'][packageVersion][dependency] = dependencyVersion.toString();
  }

  if (jsonContent.data.hasOwnProperty(name)){
    var result = {...jsonContent['data'][name]['dependencies'][packageVersion], ...data[name]['dependencies'][packageVersion]};
    if (rebuild || !jsonContent['data'][name]['dependencies'].hasOwnProperty(packageVersion)) {
      jsonContent['data'][name]['dependencies'][packageVersion] = result;
      if (description) {
        jsonContent['data'][name]['description'] = data[name]['description'];
      }
      if (repository) {
        jsonContent['data'][name]['repository'] = data[name]['repository'];
      }
      if(category) {
        jsonContent['data'][name]['category'] = data[name]['category'];
      }
    }
  } else {
    var result = {...jsonContent['data'], ...data};
    jsonContent.data = result;
  }

  fs.writeFileSync(jsonTemplateLoc, JSON.stringify(jsonContent, null, 2) + '\n');
}

export function add(args: CreateArgs) {
  const nOptions = Object.keys(args).length - 1;
  if (nOptions < 3) return false;
  if (nOptions && !Object.keys(args).slice(1).every(op =>
    ['package', 'description', 'ver', 'dep', 'depver', 'repository', 'category', 'rebuild'].includes(op))) return false;

  let repository
  if (args.repository !== undefined) {
    repository = JSON.parse(args.repository);
  } else {
    repository = undefined;
  }

  addPackageVersion(args.package, args.description, args.ver, args.dep, args.depver, repository, args.category, args.rebuild);
  return true;
}

interface CreateArgs {
  _: string[],
  package: string,
  ver: string,
  rebuild?: boolean,
  repository?: string,
  description?: string,
  dep?: string,
  depver?: string,
  category?: string
}
