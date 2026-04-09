import styles from './styles.module.css';
import React from 'react';
import Link from '@docusaurus/Link';

export default function PositionCard({name, type, location, href}) {
  return (
    <div className="col-lg-6 col-md-6 col-sm-12 my-3">
      <Link to={href} style={{textDecoration: 'none'}}>
        <div className="card border-0 shadow-sm h-100" style={{transition: 'box-shadow 0.15s ease', cursor: 'pointer'}}
          onMouseEnter={e => e.currentTarget.style.boxShadow = '0 4px 16px rgba(0,0,0,0.12)'}
          onMouseLeave={e => e.currentTarget.style.boxShadow = ''}>
          <div className="card-body text-left pt-0">
            <div className="mt-3">
              <div className="h5 mt-3 pb-2 text-dark">{name}</div>
              <div className="d-flex flex-row">
                <span className={'badge badge-pill mr-2 ' + styles.badgedg}>{type}</span>
                <span className={'badge badge-pill badge-info ' + styles.badgedg}>{location}</span>
              </div>
            </div>
          </div>
        </div>
      </Link>
    </div>
  );
}
