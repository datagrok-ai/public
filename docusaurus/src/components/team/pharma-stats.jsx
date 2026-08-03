import React from 'react';
import styles from './styles.module.css';

const stats = [
  { value: '100+', label: 'combined years in\npharma R&D' },
  { value: '7',    label: 'PhDs\non the team' },
  { value: '80+',  label: 'scientific papers\npublished' },
  { value: '10',   label: 'patented\nsoftware programs' },
];

export default function PharmaStats() {
  return (
    <section style={{background: 'var(--ifm-color-emphasis-100, #f8f9fa)', borderTop: '1px solid #e9ecef', borderBottom: '1px solid #e9ecef'}} className="py-4 mb-2">
      <div className="container">
        <h6 className="text-center text-uppercase mb-4" style={{letterSpacing: '0.12em', color: 'var(--grey-4, #6c757d)', fontSize: '0.75rem'}}>
          Why Pharma Trusts Us
        </h6>
        <div className="row justify-content-center mb-3">
          {stats.map(({value, label}) => (
            <div key={value} className="col-6 col-sm-3 text-center mb-3 mb-sm-0">
              <div style={{fontSize: '1.9rem', fontWeight: 700, lineHeight: 1, color: 'var(--ifm-color-primary-darkest, #205d3b)'}}>
                {value}
              </div>
              <div style={{fontSize: '0.78rem', color: 'var(--grey-4, #6c757d)', marginTop: '0.3rem', whiteSpace: 'pre-line', lineHeight: 1.35}}>
                {label}
              </div>
            </div>
          ))}
        </div>
        <div className="text-center mt-4">
          <a href="https://datagrok.ai/company/careers" className="btn btn-primary px-4">JOIN OUR TEAM</a>
        </div>
      </div>
    </section>
  );
}
