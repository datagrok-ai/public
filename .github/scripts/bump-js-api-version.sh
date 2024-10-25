#!/bin/bash

# Function to bump the version
bump_version() {
  local version=$1
  local level=$2

  local major=$(echo "$version" | cut -d. -f1)
  local minor=$(echo "$version" | cut -d. -f2)
  local patch=$(echo "$version" | cut -d. -f3)

  case "$level" in
    major)
      ((major++))
      minor=0
      patch=0
      ;;
    minor)
      ((minor++))
      patch=0
      ;;
    patch)
      ((patch++))
      ;;
    *)
      echo "Invalid bump level: $level"
      exit 1
      ;;
  esac

  echo "$major.$minor.$patch"
}

# Check for arguments
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <new-datagrok-version> <version-bump-level>"
  exit 1
fi

new_datagrok_version=$1
bump_level=$2

bump_json_version() {
  local package_json=$1

  # Bump package version
  current_version=$(jq -r '.version' "$package_json")
  name=$(jq -r '.name' "$package_json")
  npm_version="$(curl --retry 3 -s "https://registry.npmjs.org/${name}/${current_version}" | jq -r '.? | select( has("version") == true ).version')"
  if [[ $npm_version == "${current_version}" ]] || grep -q '\-rc' <<<"${current_version}"; then
      new_version=$(bump_version "$current_version" "$bump_level")
      jq --arg version "$new_version" '.version = $version' "$package_json" > temp.json && mv temp.json "$package_json"
  fi
}

# Find all packages with "-rc" or local relative path in dependencies, excluding node_modules
replace_datagrok_api() {
    local package_json=$1

    unpublished_jsapi=$(jq '(. | select( has("dependencies") == true ).dependencies) * (. | select( has("devDependencies") == true ).devDependencies) | select( has("datagrok-api") == true )."datagrok-api" | select(test("\\.\\./|.*-rc") == true)' "$package_json")
    if [ -n "${unpublished_jsapi}" ]; then
#      matched_packages+=("$package_json")
      # Replace datagrok-api dependency if it's using a relative path or an -rc version
      jq --arg version "^$new_datagrok_version" '.dependencies["datagrok-api"] = $version' "$package_json" > temp.json && mv temp.json "$package_json"

      bump_json_version "$package_json"

      echo "$package_json"
    fi
}
matched_packages=()
for package_json in $(find libraries packages -name 'package.json' -not -path '*/node_modules/*'); do
    matched="$(replace_datagrok_api "$package_json")"
    if [ -n "${matched}" ]; then
        matched_packages+=("${matched}")
        echo "Updated $matched to new datagrok-api and version $(jq -r '.version' "$matched")"
    fi
done

# Update dependencies where matched packages are used, excluding node_modules and preventing updates to the package's own dependency
replace_others() {
    local package_json=$1
    shift
    local matches=("$@")

    package_self_name=$(jq -r '.name' "$package_json")

    for matched_package in "${matches[@]}"; do

      package_name=$(jq -r '.name' "$matched_package")
      package_version=$(jq -r '.version' "$matched_package")

      # Only update if the dependency already exists in the package.json and avoid self-update
      levels=("dependencies" "devDependencies")
      for level in "${levels[@]}"; do
          if [[ "$package_self_name" != "$package_name" ]] && [ -n "$(jq --arg name "$package_name" --arg level "$level" '. | select( has($level) == true ) | .[$level] | select( has($name) == true )' "$package_json")" ]; then
            jq --arg name "$package_name" --arg version "^$package_version" --arg level "$level" '.[$level][$name] = $version' "$package_json" > temp.json && mv temp.json "$package_json"
            bump_json_version "$package_json"
            echo "$package_json"
          fi
      done
    done
}
for package_json in $(find libraries packages -name 'package.json' -not -path '*/node_modules/*'); do
  replace_others "$package_json" "${matched_packages[@]}"
done
