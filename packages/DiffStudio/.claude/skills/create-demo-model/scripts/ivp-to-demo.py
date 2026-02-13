"""Full-cycle exporter: IVP file -> demo TS + package.ts model + test.

Usage:
    python ivp-to-demo.py <ivp_file_path> [icon_png_path]

Examples:
    python ivp-to-demo.py files/library/pk-pd.ivp files/icons/pk-pd.png
    python ivp-to-demo.py files/library/pollution.ivp

Steps performed:
  1. Creates src/demo/<model-name>.ts with the model definition
  2. Adds import and @grok.decorators.model method to src/package.ts
  3. Adds import and testTemplate call to src/tests/demo-models-tests.ts
"""

import re
import sys
import os


def parse_ivp_header(ivp_text: str) -> tuple[str, str | None]:
    """Extract model name and description from IVP text.

    Returns:
        (name, description) where description may be None.
    """
    name = None
    description = None

    for line in ivp_text.splitlines():
        stripped = line.strip()

        if stripped.startswith('#name:'):
            name = stripped[len('#name:'):].strip()
        elif stripped.startswith('#description:'):
            description = stripped[len('#description:'):].strip()
        elif stripped.startswith('#') and name is not None:
            if not stripped.startswith('#description:'):
                break

    if name is None:
        raise ValueError('IVP file must contain a #name field')

    return name, description


def name_to_kebab(name: str) -> str:
    """Convert model name to kebab-case for file naming.

    'Acid Production' -> 'acid-production'
    'PK-PD' -> 'pk-pd'
    'Pollution' -> 'pollution'
    """
    return re.sub(r'[\s_]+', '-', name).lower()


def name_to_upper_snake(name: str) -> str:
    """Convert model name to UPPER_SNAKE_CASE for constant naming.

    'Acid Production' -> 'ACID_PRODUCTION'
    'PK-PD' -> 'PK_PD'
    'Pollution' -> 'POLLUTION'
    """
    return re.sub(r'[\s-]+', '_', name).upper()


def name_to_camel(name: str) -> str:
    """Convert model name to camelCase for method naming.

    'Acid Production' -> 'acidProduction'
    'PK-PD' -> 'pkPd'
    'Pollution' -> 'pollution'
    """
    parts = re.split(r'[\s-]+', name)
    return parts[0].lower() + ''.join(p.capitalize() for p in parts[1:])


def generate_demo_ts(ivp_content: str, name: str, description: str | None) -> str:
    """Generate the TypeScript demo file content."""
    upper_name = name_to_upper_snake(name)
    const_name = f'{upper_name}_MODEL'
    info_export = f'{upper_name}_MODEL_INFO'

    # Build the INFO model description line
    info_model_text = description if description else name

    return f"""import {{LINK}} from '../ui-constants';
import {{ModelInfo}} from '../model';

export const {const_name} = `{ivp_content}`;

const UI_OPTS = {{
  inputsTabDockRatio: 0.17,
  graphsDockRatio: 0.85,
}};

const INFO = `# Model
{info_model_text}
# Try
Interactive results based on input changes.
# Performance
Nonlinear systems of differential equations are solved within milliseconds.
# No-code
[Diff Studio](${{LINK.DIF_STUDIO}})
enables the creation of complex models without writing code.
# Learn more
* [Sensitivity analysis](${{LINK.SENS_AN}})
* [Parameter optimization](${{LINK.FITTING}})`;

export const {info_export}: ModelInfo = {{
  equations: {const_name},
  uiOptions: UI_OPTS,
  info: INFO,
}};
"""


def resolve_path(path: str, project_root: str) -> str | None:
    """Resolve a file path: try as-is (relative to cwd), then relative to project root.

    Returns absolute path if found, None otherwise.
    """
    # Try relative to cwd
    abs_path = os.path.abspath(path)
    if os.path.isfile(abs_path):
        return abs_path

    # Try relative to project root
    abs_path = os.path.join(project_root, path)
    if os.path.isfile(abs_path):
        return abs_path

    return None


def update_package_ts(package_path: str, name: str, description: str | None, icon_rel: str | None):
    """Add import and model method to package.ts."""
    kebab = name_to_kebab(name)
    upper = name_to_upper_snake(name)
    camel = name_to_camel(name)
    info_name = f'{upper}_MODEL_INFO'
    desc = (description if description else name).replace("'", "\\'")

    with open(package_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # Check for duplicates
    if info_name in content:
        print(f'  Skipped package.ts: {info_name} already present')
        return

    # --- Insert import ---
    import_pattern = re.compile(r"^import \{.*_MODEL_INFO\} from '\./demo/.*';$", re.MULTILINE)
    matches = list(import_pattern.finditer(content))
    if not matches:
        raise ValueError('Cannot find MODEL_INFO import lines in package.ts')

    insert_pos = matches[-1].end()
    import_line = f"\nimport {{{info_name}}} from './demo/{kebab}';"
    content = content[:insert_pos] + import_line + content[insert_pos:]

    # --- Insert model method ---
    # Insert before the runModel function's decorator
    marker = "  @grok.decorators.func({\n    description: 'Run model with Diff Studio UI',"
    marker_pos = content.find(marker)
    if marker_pos < 0:
        raise ValueError("Cannot find 'Run model with Diff Studio UI' marker in package.ts")

    # Build the model method block
    icon_line = f"\n    icon: '{icon_rel}'," if icon_rel else ''
    method_block = f"""  @grok.decorators.model({{
    name: '{name}',
    description: '{desc}',{icon_line}
  }})
  static async {camel}(): Promise<void> {{
    const model = new Model({info_name});
    await model.run();
  }}

"""

    content = content[:marker_pos] + method_block + content[marker_pos:]

    with open(package_path, 'w', encoding='utf-8', newline='\n') as f:
        f.write(content)

    print(f'  Updated: {package_path}')


def update_tests(tests_path: str, name: str):
    """Add import and test call to demo-models-tests.ts."""
    kebab = name_to_kebab(name)
    upper = name_to_upper_snake(name)
    info_name = f'{upper}_MODEL_INFO'

    with open(tests_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # Check for duplicates
    if info_name in content:
        print(f'  Skipped tests: {info_name} already present')
        return

    # --- Insert import ---
    import_pattern = re.compile(r"^import \{.*_MODEL_INFO\} from '\.\./demo/.*';$", re.MULTILINE)
    matches = list(import_pattern.finditer(content))
    if not matches:
        raise ValueError('Cannot find MODEL_INFO import lines in demo-models-tests.ts')

    insert_pos = matches[-1].end()
    import_line = f"\nimport {{{info_name}}} from '../demo/{kebab}';"
    content = content[:insert_pos] + import_line + content[insert_pos:]

    # --- Insert test call ---
    marker = '}); // Demo models'
    marker_pos = content.find(marker)
    if marker_pos < 0:
        raise ValueError("Cannot find '}); // Demo models' marker in tests file")

    test_line = f"  testTemplate('{kebab}', {info_name}.equations);\n"
    content = content[:marker_pos] + test_line + content[marker_pos:]

    with open(tests_path, 'w', encoding='utf-8', newline='\n') as f:
        f.write(content)

    print(f'  Updated: {tests_path}')


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    ivp_path = sys.argv[1]
    icon_arg = sys.argv[2] if len(sys.argv) > 2 else None

    # Script is at .claude/skills/create-demo-model/scripts/ivp-to-demo.py
    # Project root is 4 levels up
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, '..', '..', '..', '..'))

    # Resolve IVP file
    ivp_abs = resolve_path(ivp_path, project_root)
    if ivp_abs is None:
        print(f'Error: file not found: {ivp_path}')
        sys.exit(1)

    # Read IVP content
    with open(ivp_abs, 'r', encoding='utf-8') as f:
        ivp_content = f.read().rstrip()

    name, description = parse_ivp_header(ivp_content)
    kebab_name = name_to_kebab(name)

    # Resolve icon path
    icon_rel = None
    if icon_arg:
        icon_abs = resolve_path(icon_arg, project_root)
        if icon_abs is not None:
            icon_rel = os.path.relpath(icon_abs, project_root).replace('\\', '/')
        else:
            print(f'  Warning: icon file not found: {icon_arg}, skipping icon')

    # Step 1: Create demo TS file
    demo_path = os.path.join(project_root, 'src', 'demo', f'{kebab_name}.ts')
    os.makedirs(os.path.dirname(demo_path), exist_ok=True)
    ts_content = generate_demo_ts(ivp_content, name, description)
    with open(demo_path, 'w', encoding='utf-8', newline='\n') as f:
        f.write(ts_content)
    print(f'  Created: {demo_path}')

    # Step 2: Update package.ts
    package_path = os.path.join(project_root, 'src', 'package.ts')
    update_package_ts(package_path, name, description, icon_rel)

    # Step 3: Update tests
    tests_path = os.path.join(project_root, 'src', 'tests', 'demo-models-tests.ts')
    update_tests(tests_path, name)

    print(f'\nDone! Model "{name}" fully registered.')


if __name__ == '__main__':
    main()
