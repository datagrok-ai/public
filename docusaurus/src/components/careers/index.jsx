import styles from './styles.module.css';
import React from 'react';

export function LayoutHeader ({meta}){
    return(
        <div className='row'>
                <div className="col-lg-12 mb-12 px-5 pb-3">
                    <div className='display-5 text-secondary'>{meta.title}</div>
                    <h4 className='mt-2 pb-2'>{meta.description}</h4>
                </div>
        </div>
        
    );
}
