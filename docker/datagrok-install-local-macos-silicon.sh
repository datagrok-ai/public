#!/bin/bash

# This script automates the Datagrok local installation and running
# To see additional actions, run "./datagrok-install-local help"

set -o errexit

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
RESET='\033[0m'
timeout=30

# Any tag of datagrok/datagrok: `latest` for the newest stable release, an exact `X.Y.Z`,
# `bleeding-edge` to track master, or an `X.Y.Z-rc` candidate.
datagrok_default_version="1.27.9"

# A compose file and the images it starts must come from the same release: the compose file
# carries the GROK_PARAMETERS that the datlas build inside datagrok/datagrok parses, and the
# two evolve together. So the compose file is not taken from whatever master happens to hold
# — it is taken from the branch the requested image was actually built from, which the image
# records in its `branch` label. That works for every tag, including moving ones like
# `latest`, and it means the pairing is derived from the image rather than guessed from the
# tag's name.
DATAGROK_VERSION="${DATAGROK_VERSION:-$datagrok_default_version}"

compose_config_name="localhost.macos-silicon.docker-compose.yaml"
datagrok_local_url="http://localhost:8080/"
datagrok_image="datagrok/datagrok:${DATAGROK_VERSION}"

script_name="$0"
script_dir=$(cd -- "$(dirname -- "$0")" && pwd)
compose_config_path="${script_dir}/${compose_config_name}"

# Filled in by resolve_release: the git ref to take the compose file from, the release that
# ref's compose file declares, and whether the copy on disk is older than the tag it carries.
resolved_ref=""
resolved_release=""
image_stale=false

function message() {
    echo -e "${YELLOW}$1${RESET}"
}

function error() {
    echo -e "${RED}!!!$1!!!${RESET}"
}

function clean_label() {
    local v
    v=$(tr -d '\r' | tail -1)
    [ "$v" = "<no value>" ] && v=""
    echo "$v"
}

# What the registry says, which for a moving tag like `latest` is the only authority — a copy
# on disk may point at an older release. imagetools reads the image config straight from the
# registry, so this costs no download.
function registry_branch() {
    docker buildx imagetools inspect "${datagrok_image}" \
        --format '{{ index .Image.Config.Labels "branch" }}' 2>/dev/null | clean_label
}

function local_branch() {
    docker image inspect "${datagrok_image}" \
        --format '{{index .Config.Labels "branch"}}' 2>/dev/null | clean_label
}

# Work out which branch built the image this tag names. Registry first; a local image is the
# fallback for when the registry can't be reached, so an existing install still starts offline.
function resolve_release() {
    [ -n "${resolved_ref}" ] && return 0
    local branch local_b
    branch=$(registry_branch)
    local_b=$(local_branch)
    # A local copy older than the tag now points at has to be replaced: running it against the
    # newer release's compose file is exactly the pairing this all exists to prevent.
    if [ -n "$branch" ] && [ -n "$local_b" ] && [ "$branch" != "$local_b" ]; then
        image_stale=true
    fi
    if [ -z "$branch" ]; then
        branch="$local_b"
    fi
    if [ -z "$branch" ]; then
        message "Resolving ${DATAGROK_VERSION}"
        docker pull "${datagrok_image}" >/dev/null || {
            error "Could not pull ${datagrok_image}"
            message "Check that '${DATAGROK_VERSION}' is a published tag of datagrok/datagrok."
            exit 255
        }
        branch=$(local_branch)
    fi
    case "$branch" in
        release/*)
            resolved_ref="$branch"
            resolved_release="${branch#release/}"
            ;;
        "")
            # Images built before the label existed. Fall back to reading the tag, which is
            # right for the two shapes those old tags came in.
            if [ "${DATAGROK_VERSION}" = "bleeding-edge" ]; then
                resolved_ref="master"; resolved_release="bleeding-edge"
            elif [[ ${DATAGROK_VERSION} =~ ^[0-9]+\.[0-9]+\.[0-9]+ ]]; then
                resolved_ref="release/${DATAGROK_VERSION%%-*}"
                resolved_release="${DATAGROK_VERSION%%-*}"
            else
                error "Cannot tell which release ${datagrok_image} belongs to"
                message "The image carries no 'branch' label and '${DATAGROK_VERSION}' is not a version."
                message "Install an explicit version instead: DATAGROK_VERSION=1.27.7"
                exit 255
            fi
            message "Note: ${datagrok_image} has no branch label; assuming ${resolved_release}."
            ;;
        *)
            # A build from master or a feature branch: those track master's compose file.
            resolved_ref="master"
            resolved_release="bleeding-edge"
            ;;
    esac
    if [ "${DATAGROK_VERSION}" != "${resolved_release}" ]; then
        message "${DATAGROK_VERSION} is ${resolved_release} (built from ${resolved_ref})"
    fi
}

function compose_url() {
    echo "https://raw.githubusercontent.com/datagrok-ai/public/${resolved_ref}/docker/${compose_config_name}"
}

# The release a compose file on disk was cut from, or empty if it predates the marker.
function compose_release() {
    [ -f "${compose_config_path}" ] || return 0
    grep -E '^x-datagrok-release:' "${compose_config_path}" | head -1 | sed 's/^[^:]*:[[:space:]]*//' | tr -d '\r'
}

# Refuse to start a compose file from one release against images from another. The comparison
# is against the release the image resolved to, not the tag the user typed, so a moving tag
# like `latest` is checked against what it currently points at. Files cut before the marker
# existed have nothing to compare; those are only warned about, since a compose file taken
# from release/X.Y.Z already matches X.Y.Z by construction.
function verify_compose_release() {
    local found
    found=$(compose_release)
    if [ -z "$found" ]; then
        message "Note: ${compose_config_name} predates release binding; assuming ${resolved_release}."
        return 0
    fi
    if [ "$found" != "${resolved_release}" ]; then
        error "Version mismatch"
        message "${compose_config_name} is built for '${found}', but ${DATAGROK_VERSION} is '${resolved_release}'."
        message "Running them together makes the Datagrok server fail to start."
        message "Delete ${compose_config_path} and re-run this script to download the matching file."
        exit 255
    fi
}

function download_compose() {
    message "Downloading Datagrok config file for ${resolved_release}"
    curl -fsSL -o "${compose_config_path}" "$(compose_url)" || {
        error "Could not download the config file"
        message "Tried: $(compose_url)"
        message "${datagrok_image} was built from ${resolved_ref}, but that branch has no compose file."
        exit 255
    }
    verify_compose_release
}

function user_query_yn() {
    echo -ne "${YELLOW}"
    read -r -p "${1} (y/N)" answer
    echo -ne "${RESET}"
    if [[ $answer =~ ^(Y|y) ]]; then
        return 0 # 0 = True
    else
        return 1 # 1 = False
    fi
}

function count_down() {
    echo -n "Waiting:"
    for ((i = ${1}; i > 0; i--)); do
        echo -n " $i"
        sleep 1
    done
    echo " 0"
}

function check_docker() {
    if [ ! -x "$(command -v docker)" ]; then
        error "Docker engine is not installed"
        message "Please install Docker and Docker Compose plugin by manual: https://docs.docker.com/engine/install/"
        exit 255
    fi
}

function check_docker_daemon() {
    docker info >/dev/null || {
        error "Docker daemon is not running"
        message "Please launch Docker Desktop application"
        exit 255
    }
}

function check_installation() {
    if [ ! -f "${compose_config_path}" ]; then
        return 1 # False
    else
        docker_image=$(docker images -q "datagrok/datagrok")
        if [ ! -n "$docker_image" ]; then
            return 1 # False
        fi
    fi
    return 0 # True
}

function datagrok_install() {
    resolve_release
    # A file left over from another release is replaced, not reused.
    if [ ! -f "${compose_config_path}" ] || [ "$(compose_release)" != "${resolved_release}" ]; then
        download_compose
    fi
    verify_compose_release
    message "Pulling Datagrok ${resolved_release} images (this can take a while depending on your Internet connection speed)"
    docker compose -f "${compose_config_path}" --profile all pull
}

function datagrok_start() {
    if ! check_installation; then
        datagrok_install
    fi
    message "Starting Datagrok containers"
    update_installation=false
    if [ ! -z "$1" ] && [ "$1"=="update" ]; then
        update_installation=true
    fi
    # Checking do we need tu run update
    if [ "$update_installation" = true ] ; then
        resolve_release
        download_compose
        message "Updating Datagrok to ${resolved_release}"
        docker compose -f "${compose_config_path}" --project-name datagrok --profile all pull
        docker compose -f "${compose_config_path}" --project-name datagrok --profile all  up -d --force-recreate
    else
        resolve_release
        verify_compose_release
        # The tag moved since this image was pulled; take the new one so the running server
        # matches the compose file about to configure it.
        if [ "${image_stale}" = true ]; then
            message "${DATAGROK_VERSION} has moved to ${resolved_release}; pulling"
            docker compose -f "${compose_config_path}" --profile all pull
        fi
        docker compose -f "${compose_config_path}" --project-name datagrok --profile all up -d
    fi
    message "Waiting while the Datagrok server is starting"
    echo "When the browser opens, use the following credentials to log in:"
    echo "------------------------------"
    echo -ne "${GREEN}"
    echo "Login:    admin"
    echo "Password: admin"
    echo -ne "${RESET}"
    echo "------------------------------"
    echo "If you see the message 'Datagrok server is unavaliable' just wait for a while and reload the web page "
    count_down ${timeout}
    message "Running browser"
    open ${datagrok_local_url}
    message "If the browser hasn't open, use the following link: $datagrok_local_url"
    message "To extend Datagrok fucntionality, install extension packages on the 'Manage -> Packages' page"
    if [ "$update_installation" = true ] ; then
        message "Removing old images"
        docker image prune -f
    fi
}

function datagrok_stop() {
    if ! check_installation; then
        message "The Datagrok installation not found. Nothing to stop."
        exit 255
    fi
    message "Stopping Datagrok containers"
    docker compose -f "${compose_config_path}" --project-name datagrok --profile all stop
    docker rm -f $(docker ps --format {{.Names}} | grep datagrok)
}

function datagrok_reset() {
    if ! check_installation; then
        message "The Datagrok installation not found. Nothing to reset."
        exit 255
    fi
    if user_query_yn "This action will stop Datagrok, remove all user settings, data and installed packages. Are you sure?"; then
        docker compose -f "${compose_config_path}" --project-name datagrok --profile all stop
        docker compose -f "${compose_config_path}" --project-name datagrok --profile all down --volumes
    fi
}

function datagrok_purge() {
    if ! check_installation; then
        message "The Datagrok installation not found. Nothing to remove."
        exit 255
    fi
    if user_query_yn "This action will stop Datagrok and COMPLETELY remove Datagrok installation. Are you sure?"; then
        docker compose -f "${compose_config_path}" --project-name datagrok --profile all down --volumes
        docker rmi $(docker images -q datagrok/*)
    fi
}

# === Main part of the script starts from here ===
check_docker
check_docker_daemon

case "$1" in
install) datagrok_install ;;
start) datagrok_start ;;
stop) datagrok_stop ;;
reset) datagrok_reset ;;
update) datagrok_start update ;;
purge) datagrok_purge ;;
version)
    resolve_release
    echo "${DATAGROK_VERSION} -> ${resolved_release} (${resolved_ref})"
    ;;
help | "-h" | "--help")
    echo "usage: $script_name install|start|stop|update|reset|purge|version" >&2
    echo "" >&2
    echo "  DATAGROK_VERSION is any tag of datagrok/datagrok (default ${datagrok_default_version}):" >&2
    echo "    latest         newest stable release" >&2
    echo "    <X.Y.Z>        that exact release" >&2
    echo "    <X.Y.Z>-rc     that release candidate" >&2
    echo "    bleeding-edge  latest build from master" >&2
    echo "" >&2
    echo "  The compose file is taken from the branch the chosen image was built from," >&2
    echo "  so the configuration always matches the images." >&2
    exit 1
    ;;
*) datagrok_start ;;
esac
