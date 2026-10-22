import json
from pathlib import Path
from typing import Any

from ruamel.yaml import YAML

from .protocol import Settings

HERE = Path(__file__).parent


def load_yaml(name: str) -> Any:
    return YAML(typ='safe').load(HERE / name)


def write_home(home: Path, settings: Settings) -> None:
    config = load_yaml('hermes.yaml')
    config['skills']['external_dirs'] = [settings['skillsDir']]
    config['compression']['threshold_tokens'] = settings['compressionThresholdTokens']
    # As JSON, which is YAML, a string such as "off" reaches Hermes as a string whatever YAML version it reads.
    (home / 'config.yaml').write_text(json.dumps(config))
    plugins = home / 'plugins'
    plugins.mkdir(mode=0o700, parents=True, exist_ok=True)
    link = plugins / 'grokky'
    link.unlink(missing_ok=True)
    link.symlink_to(HERE / 'hermes_plugin')
