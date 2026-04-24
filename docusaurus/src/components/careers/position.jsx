import styles from './styles.module.css';
import React from 'react';

export default function profileCard({meta}) {
    return (
        <div id={meta.id} className={'col-lg-6 col-md-6 col-sm-12 my-3'}>
            <div className='card hover-shadow border-0 border-light shadow-sm h-100'>
                <a href={'/company/careers/' + meta.id} className='btn link'>
                    <div className='card-body text-left pt-0'>
                        <div className='mt-3'>
                            <div className='h5 mt-3 pb-2'>{meta.name}</div>
                            <div className='d-flex flex-row'>
                                <span className={'badge badge-pill mr-2 ' + styles.badgedg}>{meta.type}</span>
                                <span className={'badge badge-pill badge-info ' + styles.badgedg}>{meta.location}</span>
                            </div>
                        </div>
                    </div>
                </a>
            </div>
        </div>
    );
}
