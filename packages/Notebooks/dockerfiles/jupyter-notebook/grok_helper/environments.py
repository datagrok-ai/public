#!/usr/bin/env python
import os
import re
import subprocess
import time
import yaml

from utilites import call_process

BASE_DEPENDENCIES = ['ipykernel', 'matplotlib', 'pandas', 'pip', 'requests', 'snappy<1.2', 'pyarrow']
PIP_DEPENDENCIES = ['datagrok-api']
BASE_DEPENDENCIES_SET = frozenset(BASE_DEPENDENCIES)
CONDA_HOME = os.environ.get('CONDA_HOME', '/home/grok/conda')
TEMPLATE_ENV_PREFIX = '_base_py'

R_BASE_DEPENDENCIES = [
    'r-base', 'r-irkernel', 'r-renv', 'r-httr', 'r-repr', 'r-data.table',
    'r-r.utils', 'r-fs', 'r-arrow', 'r-jsonlite',
]

LANGUAGE_CONFIG = {
    'python': {
        'base_dependencies': BASE_DEPENDENCIES,
        'pip_dependencies': PIP_DEPENDENCIES,
        'extra_channels': [],
    },
    'r': {
        'base_dependencies': R_BASE_DEPENDENCIES,
        'pip_dependencies': [],
        'extra_channels': ['bioconda'],
    },
}

def _extract_pip_from_dependencies(dependencies):
    conda_deps = []
    pip_packages = []
    for dep in dependencies:
        if isinstance(dep, dict) and 'pip' in dep:
            pip_packages.extend(dep['pip'])
        else:
            conda_deps.append(dep)
    return conda_deps, pip_packages

def _validate_dependencies(dependencies):
    for dep in dependencies:
        if isinstance(dep, dict):
            for key in dep:
                if key != 'pip':
                    raise ValueError(
                        f"Unsupported dependency type '{key}:' in environment YAML. "
                        f"Use 'pip:' for PyPI packages."
                    )

class Environments(object):
    # Path to cache directory
    path = 'environments'

    def __init__(self, path='environments'):
        if not os.path.exists(path):
            os.mkdir(path)
        self.path = path

    def list(self) -> list:
        """
        Returns list existing environments.

        :return: List existing environments.
        """
        log = call_process(['micromamba', 'env', 'list'], self.path).split('\n')
        del log[0:1]
        envs = []
        for line in log:
            if not line.startswith('#') and len(line.strip()) > 0:
                match = re.search('^ *([\w\-]+)', line)
                if match is not None:
                    envs.append(match.group(1))
        return envs


    def _get_python_version(self, data):
        for dep in (data.get('dependencies') or []):
            if isinstance(dep, str) and dep.startswith('python='):
                parts = dep.split('=')[1].split('.')
                if len(parts) >= 2:
                    return (int(parts[0]), int(parts[1]))
        return None

    def _get_template_tar_path(self, data):
        for dep in (data.get('dependencies') or []):
            if isinstance(dep, str) and dep.startswith('python='):
                ver = dep.split('=')[1]
                major_minor = '.'.join(ver.split('.')[:2])
                tar_path = os.path.join(CONDA_HOME, 'envs', TEMPLATE_ENV_PREFIX + major_minor + '.tar')
                if os.path.isfile(tar_path):
                    return tar_path
        return None

    def _normalize_environment(self, environment, language):
        config = LANGUAGE_CONFIG.get(language, LANGUAGE_CONFIG['python'])
        data = yaml.safe_load(environment)
        if data['dependencies'] is None:
            data['dependencies'] = []
        data['dependencies'].extend(config['base_dependencies'])
        _validate_dependencies(data['dependencies'])
        conda_deps, yaml_pip = _extract_pip_from_dependencies(data['dependencies'])
        data['dependencies'] = conda_deps
        python_version = self._get_python_version(data)
        skip_pip = language == 'python' and python_version is not None and python_version < (3, 8)
        data['pip_dependencies'] = ([] if skip_pip else list(config['pip_dependencies'])) + yaml_pip
        return data, config

    def exists(self, id: str, environment: str, language: str = 'python') -> bool:
        result = subprocess.run(['jupyter', 'kernelspec', 'list'], capture_output=True, text=True)
        if id not in result.stdout:
            return False
        environment_path = os.path.join(self.path, id + '-extended.yaml')
        if os.path.isfile(environment_path):
            with open(environment_path, "r") as f:
                current_env = yaml.safe_load(f.read())
                check_env, _ = self._normalize_environment(environment, language)
                current_env.setdefault('pip_dependencies', [])
                return current_env['dependencies'] == check_env['dependencies'] and \
                    current_env.get('pip_dependencies') == check_env.get('pip_dependencies')
        return False


    def create(self, id: str, environment: str, language: str = 'python'):
        """
        Creates environment.

        :param id: Environment ID (name).
        :param environment: Environment configuration.
        :param language: Script language ('python' or 'r').
        :return: Log.
        """
        env_file = os.path.join(self.path, id + '.yaml')
        file = open(env_file, 'w+')
        file.write(environment)
        file.close()

        data, config = self._normalize_environment(environment, language)

        channels = data.get('channels', [])
        if 'conda-forge' not in channels:
            channels.append('conda-forge')
        for ch in config['extra_channels']:
            if ch not in channels:
                channels.append(ch)
        data['channels'] = channels

        extended_env_file = os.path.join(self.path, id + '-extended.yaml')
        with open(extended_env_file, 'w+') as extended_file:
            yaml.dump(data, extended_file, default_flow_style=False)

        extra_deps = [d for d in data['dependencies']
                     if d not in BASE_DEPENDENCIES_SET and not (isinstance(d, str) and d.startswith('python='))]

        template_tar = self._get_template_tar_path(data)
        env_path = os.path.join(CONDA_HOME, 'envs', id)

        script_file = os.path.join(self.path, id + '.sh')
        script = \
            'set -e -o pipefail\n' + \
            'source /etc/profile.d/micromamba.sh\n' + \
            'JUPYTER_BIN_DIR=$(dirname "$(which jupyter)")\n' + \
            'if micromamba env list -q | grep -q ' + id + '; then\n' + \
            'echo "Environment already exists, updating"\n' + \
            'micromamba install -y --name ' + id + ' --file=' + extended_env_file + '\n' + \
            'else\n'

        if template_tar:
            # Extract pre-built template tarball (sequential I/O, fast on overlay2).
            # Then fix shebangs and symlinks that reference the template env path.
            # Use grep -rIl (-I skips binary files) to avoid corrupting ELF binaries.
            template_env_dir = os.path.splitext(template_tar)[0]
            script += \
                'mkdir -p ' + env_path + '\n' + \
                'tar xf ' + template_tar + ' -C ' + env_path + '\n' + \
                'grep -rIl "' + template_env_dir + '" ' + env_path + '/bin/ ' + env_path + '/etc/ 2>/dev/null ' + \
                "| xargs sed -i 's|" + template_env_dir + '|' + env_path + "|g' 2>/dev/null || true\n" + \
                'find ' + env_path + ' -type l -lname "' + template_env_dir + '/*" | while read lnk; do\n' + \
                'target=$(readlink "$lnk")\n' + \
                'ln -sf "${target/' + template_env_dir + '/' + env_path + '}" "$lnk"\n' + \
                'done\n' + \
                'export LD_LIBRARY_PATH=' + env_path + '/lib:${LD_LIBRARY_PATH:-}\n'
            if extra_deps:
                script += \
                    'micromamba install -y --name ' + id + ' --file=' + extended_env_file + '\n'
        else:
            script += \
                'micromamba env create -y --name ' + id + ' --file=' + extended_env_file + '\n'

        script += \
            'fi\n' + \
            'micromamba activate ' + id + '\n'

        if data['pip_dependencies']:
            # Generate upper-bound constraints from conda-installed packages to prevent
            # uv/pip from upgrading them (e.g. numpy 1.x -> 2.x which breaks pyarrow ABI).
            # Use <= instead of == so uv can pick older compatible transitive deps
            # when requires-python metadata conflicts (e.g. zipp on Python 3.8).
            python_version = self._get_python_version(data)
            uv_python_flag = ' --python-version ' + '.'.join(str(x) for x in python_version) if python_version else ''
            script += 'PYTHON_PATH=$(which python)\n' + \
                'if [[ "$PYTHON_PATH" != *"' + id + '"* ]]; then\n' + \
                'echo "ERROR: micromamba activate failed - python resolves to $PYTHON_PATH"\n' + \
                'exit 1\n' + \
                'fi\n' + \
                'micromamba list -n ' + id + ' --json | python -c "import sys,json; print(chr(10).join(p[\\"name\\"]+\\"<=\\"+p[\\"version\\"] for p in json.load(sys.stdin) if not p[\\"name\\"].startswith(\\"_\\")))" > /tmp/' + id + '_constraints.txt\n' + \
                'uv pip install --system' + uv_python_flag + ' --constraint /tmp/' + id + '_constraints.txt ' + ' '.join(data['pip_dependencies']) + '\n' + \
                'rm -f /tmp/' + id + '_constraints.txt\n'

        if language == 'r':
            script += 'export PATH="$PATH:$JUPYTER_BIN_DIR"\n'
            script += "Rscript -e \"IRkernel::installspec(name='" + id + "', displayname='" + id + "', user=TRUE)\"\n"
        else:
            script += 'python -m ipykernel install --user --name ' + id + '\n'

        file = open(script_file, 'w+')
        file.write(script)
        file.close()

        log = call_process(['bash', script_file], self.path)
        self.wait_for_jupyter_kernel(id)
        return log


    def remove(self, id):
        """
        Removes environment.

        :param id: Environment ID (name).
        :return: Log.
        """
        return call_process(['micromamba', 'remove', '--name', id, '--all', '--yes'], self.path)


    def wait_for_jupyter_kernel(self, id: str, retries=20, wait_time=3):
        while retries > 0:
            result = subprocess.run(['jupyter', 'kernelspec', 'list'], capture_output=True, text=True)
            if id in result.stdout:
                return True
            time.sleep(3)
            retries -= 1
        return False
