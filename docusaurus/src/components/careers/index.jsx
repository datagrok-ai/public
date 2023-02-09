import styles from './styles.module.css';
import React from 'react';

export function LayoutHeader ({children, meta}){
    return(
        <div className="row">
            <div className="row py-2 align-items-center">
                <div className="col-lg-12 px-5">
                    <div className='display-5 text-secondary pb-3'>{meta.title}</div>
                    <p class="h6">{meta.description}</p>
                </div>
            </div>
        </div>
    );
}