import styles from './styles.module.css';
import React from 'react';
import Link from '@docusaurus/Link';

export default function PositionPage({children, meta}) {
  return (
    <div className="container py-5">
      <div className="row justify-content-center">
        <div className="col-lg-8 col-md-10 col-12">

          <Link to="/help/careers/default" style={{fontSize: '0.875rem', color: 'var(--ifm-color-primary-dark)'}}>
            ← All open positions
          </Link>

          <h1 className="mt-3 mb-2">{meta.name}</h1>

          <div className="d-flex flex-row mb-4">
            <span className={'badge badge-pill mr-2 ' + styles.badgedg}>{meta.type}</span>
            <span className={'badge badge-pill badge-info ' + styles.badgedg}>{meta.location}</span>
          </div>

          <div className="text-left">
            {children}
          </div>

          <div className="mt-4 pt-2">
            <a className="btn btn-primary px-4" href={'mailto:hr@datagrok.ai?subject=' + encodeURIComponent(meta.name)}>
              Apply for this position
            </a>
          </div>

        </div>
      </div>
    </div>
  );
}
