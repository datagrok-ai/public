#!/bin/sh

set -ex

RESOLVER_ADDR="$(cat /etc/resolv.conf | grep nameserver | head | awk '{print $2}')"
export RESOLVER_ADDR
envsubst '${RESOLVER_ADDR}' < /etc/nginx/nginx.conf.template > /etc/nginx/nginx.conf


nginx -g 'daemon off;' & exec /home/grok/entrypoint.sh \
    --main "jupyter notebook --debug --config=/home/grok/jupyter_notebook_config.py" \
    --helper "gunicorn grok_helper:app --timeout 900 --chdir ${GROK_HELPER_DIR} -b 0.0.0.0:5005 --access-logfile - --error-logfile - --workers 2 --threads=4 --worker-class=gthread" "$@"
