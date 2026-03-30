import fs from 'fs';
import path from 'path';
import * as color from './color-utils';

// Header tags recognized in Python function metadata comments (ported from ScriptParser.headerTags)
const headerTags = [
  'name', 'description', 'help-url', 'input', 'output', 'tags', 'sample', 'language',
  'endpoint', 'requiresServer', 'param-csrfmiddlewaretoken', 'returns', 'test', 'sidebar',
  'condition', 'top-menu', 'environment', 'require', 'editor-for', 'schedule', 'schedule.runAs',
  'reference', 'editor',
];

const paramRegex = new RegExp(`(#\\s*(${headerTags.join('|')}|meta\\.[^:]*):)\\s*(\\S+)\\s?(\\S+)?`);
const defRegex = /^\s*def\s+([a-zA-Z_][a-zA-Z0-9_]*)\s*\(/;

const TASKS_FILE = 'tasks.yaml';

const TEMPLATE_PIP_DOCKER = `FROM python:3.11-slim-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \\
    gcc \\
    libc-dev \\
    && apt-get clean \\
    && rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

RUN useradd -m datagrok

WORKDIR /app
COPY . /app
ENV PYTHONPATH="/app:\${PYTHONPATH}"
RUN chown -R datagrok:datagrok /app

# Install Python dependencies using uv
RUN uv pip install --system celery datagrok-celery-task pyyaml && \\
    if [ -f requirements.in ]; then uv pip install --system --no-cache -r requirements.in; \\
    elif [ -f requirements.txt ]; then uv pip install --system --no-cache -r requirements.txt; fi

USER datagrok

EXPOSE 8000

CMD celery -A \$DATAGROK_CELERY_NAME worker \\
    --loglevel=info \\
    --hostname=\$CELERY_HOSTNAME \\
    --concurrency=\$WORKERS_PROCESSES \\
    -Q \$TASK_QUEUE_NAME
`;

const TEMPLATE_CONDA_DOCKER = `FROM mambaorg/micromamba:latest as builder

USER root

ENV MAMBA_ROOT_PREFIX=/opt/conda \\
    MAMBA_DOCKERFILE_ACTIVATE=1 \\
    DEFAULT_ENV_NAME=myenv

COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

WORKDIR /tmp/app

COPY environment.yaml environment.yml* ./

RUN if [ -f environment.yaml ]; then sed -i 's/\\r\$//' environment.yaml; fi && \\
    if [ -f environment.yml ]; then sed -i 's/\\r\$//' environment.yml; fi

RUN ENV_FILE=\$( [ -f environment.yaml ] && echo environment.yaml || echo environment.yml ) && \\
    ENV_NAME=\$(grep '^name:' \$ENV_FILE | cut -d' ' -f2 || echo "myenv") && \\
    micromamba create -y -n \$ENV_NAME -f \$ENV_FILE && \\
    uv pip install --python /opt/conda/envs/\$ENV_NAME/bin/python celery datagrok-celery-task pyyaml && \\
    echo "\$ENV_NAME" > /tmp/env_name.txt

RUN micromamba clean --all --yes

FROM mambaorg/micromamba:latest

USER root

ENV MAMBA_ROOT_PREFIX=/opt/conda \\
    MAMBA_DOCKERFILE_ACTIVATE=1

# Create non-root user
RUN useradd -m -s /bin/bash grok

WORKDIR /app

COPY --from=builder /opt/conda /opt/conda
COPY --from=builder /tmp/env_name.txt /opt/env_name.txt

COPY . /app

RUN chown -R grok:grok /opt/conda /app /opt/env_name.txt
EXPOSE 8000
USER grok
SHELL ["/bin/bash", "-c"]

CMD micromamba run -n \$(cat /opt/env_name.txt) \\
    celery -A \$DATAGROK_CELERY_NAME worker \\
        --loglevel=info \\
        --hostname=\$CELERY_HOSTNAME \\
        --concurrency=\$WORKERS_PROCESSES \\
        -Q \$TASK_QUEUE_NAME
`;

const TEMPLATE_PYTHON_ENTRY = `import os
import importlib
import logging

import yaml
from celery import Celery
from datagrok_celery_task import DatagrokTask, Settings

import sys
sys.path.insert(0, os.getcwd())

settings = Settings(log_level=logging.DEBUG)
app = Celery(settings.celery_name, broker=settings.broker_url)

with open("${TASKS_FILE}") as f:
    tasks = yaml.safe_load(f)

if not isinstance(tasks, dict) or not isinstance(tasks['tasks'], list):
    raise RuntimeError('Incorrect format for tasks.yaml. It should contain tasks field.')

for task_def in tasks['tasks']:
    mod = importlib.import_module(task_def["module"])
    fn = getattr(mod, task_def["name"])
    app.task(name=task_def["name"], base=DatagrokTask)(fn)
`;

interface TaskEntry {
  name: string;
  module: string;
}

function getAllPyFiles(dir: string): string[] {
  const results: string[] = [];
  const entries = fs.readdirSync(dir, {withFileTypes: true});
  for (const entry of entries) {
    const fullPath = path.join(dir, entry.name);
    if (entry.isDirectory())
      results.push(...getAllPyFiles(fullPath));
    else if (entry.name.endsWith('.py'))
      results.push(fullPath);
  }
  return results;
}

function scanPythonFunctions(pythonDir: string, root: string): TaskEntry[] {
  const pyFiles = getAllPyFiles(pythonDir);
  const tasks: TaskEntry[] = [];

  for (const pyFile of pyFiles) {
    const content = fs.readFileSync(pyFile, 'utf-8');
    const lines = content.split('\n');
    let header: string[] = [];
    let isHeader = false;

    for (const line of lines) {
      if (line.length > 0 && paramRegex.test(line)) {
        if (!isHeader)
          isHeader = true;
        header.push(line);
      }
      else {
        if (isHeader && defRegex.test(line)) {
          const match = defRegex.exec(line);
          const name = match![1];
          let module = path.basename(pyFile, '.py');
          const fileDir = path.dirname(pyFile);
          if (fileDir !== root) {
            const relDir = path.relative(root, fileDir).replace(/[\\/]+/g, '.');
            if (relDir && relDir !== '.')
              module = relDir + '.' + module;
          }
          tasks.push({name, module});
          isHeader = false;
          header = [];
        }
        else if (isHeader) {
          isHeader = false;
          header = [];
        }
      }
    }
  }
  return tasks;
}

function copyDirContents(src: string, dest: string): void {
  const entries = fs.readdirSync(src, {withFileTypes: true});
  for (const entry of entries) {
    const srcPath = path.join(src, entry.name);
    const destPath = path.join(dest, entry.name);
    if (entry.isDirectory()) {
      fs.mkdirSync(destPath, {recursive: true});
      copyDirContents(srcPath, destPath);
    }
    else
      fs.copyFileSync(srcPath, destPath);
  }
}

function useConda(dir: string): boolean {
  return fs.existsSync(path.join(dir, 'environment.yaml')) ||
    fs.existsSync(path.join(dir, 'environment.yml'));
}

function deployFolder(packageDir: string, folderPath: string, dirName: string): boolean {
  const tasks = scanPythonFunctions(folderPath, folderPath);
  if (tasks.length === 0) {
    color.log(`No annotated Python functions found in ${dirName}`);
    return false;
  }

  const dockerfilesDir = path.join(packageDir, 'dockerfiles', dirName);
  fs.mkdirSync(dockerfilesDir, {recursive: true});

  // Generate Dockerfile
  const dockerfile = useConda(folderPath) ? TEMPLATE_CONDA_DOCKER : TEMPLATE_PIP_DOCKER;
  fs.writeFileSync(path.join(dockerfilesDir, 'Dockerfile'), dockerfile);

  // Generate tasks.yaml (JSON-encoded, matching server behavior)
  fs.writeFileSync(path.join(dockerfilesDir, TASKS_FILE), JSON.stringify({tasks}, null, 2));

  // Generate Celery entry point
  const celeryName = dirName.replace(/-/g, '_');
  fs.writeFileSync(path.join(dockerfilesDir, celeryName + '.py'), TEMPLATE_PYTHON_ENTRY);

  // Copy Python source files
  copyDirContents(folderPath, dockerfilesDir);

  // Copy container.json if present
  const containerJsonSrc = path.join(folderPath, 'container.json');
  const containerJsonDest = path.join(dockerfilesDir, 'container.json');
  if (fs.existsSync(containerJsonSrc) && !fs.existsSync(containerJsonDest))
    fs.copyFileSync(containerJsonSrc, containerJsonDest);

  color.log(`Generated Celery Docker artifacts in dockerfiles/${dirName}/`);
  return true;
}

export function generateCeleryArtifacts(packageDir: string): boolean {
  const pythonDir = path.join(packageDir, 'python');
  if (!fs.existsSync(pythonDir))
    return false;

  const entries = fs.readdirSync(pythonDir, {withFileTypes: true});
  const isNested = entries.length > 0 && entries.every((e) => e.isDirectory);
  let generated = false;

  if (isNested) {
    for (const entry of entries) {
      if (!entry.isDirectory())
        continue;
      const folderPath = path.join(pythonDir, entry.name);
      if (deployFolder(packageDir, folderPath, entry.name))
        generated = true;
    }
  }
  else {
    const dirName = path.basename(packageDir).toLowerCase();
    if (deployFolder(packageDir, pythonDir, dirName))
      generated = true;
  }

  return generated;
}
