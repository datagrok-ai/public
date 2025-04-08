FROM datagrok/jupyter_notebook:latest

USER root
RUN apt-get update && apt-get install -y nginx gettext && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /etc/nginx/sites-enabled /etc/nginx/conf.d

# Copy nginx configuration files
COPY nginx.conf.template /etc/nginx/nginx.conf.template
COPY grok_cvm.conf.template /etc/nginx/sites-enabled/grok_cvm.conf
COPY entrypoint.sh /etc/nginx/entrypoint.sh
COPY jupyter_notebook_config.py /home/grok/jupyter_notebook_config.py
COPY grok_helper_configuration.yaml ${GROK_HELPER_DIR}/grok_helper_configuration.yaml
RUN sed -i 's/\r$//' /etc/nginx/entrypoint.sh && chmod +x /etc/nginx/entrypoint.sh

RUN chown -R grok:grok /etc/nginx /var/lib/nginx /var/log/nginx

EXPOSE 8090

USER grok
WORKDIR '/etc/nginx'
ENTRYPOINT ["/bin/sh", "/etc/nginx/entrypoint.sh"]
