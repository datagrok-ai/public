import styles from './styles.module.css';
import React from 'react';

export function LayoutHeader ({children, meta}){
    return null;
}
export function LayoutBody ({children, meta}){
    return(

        <div className="row justify-content-md-center text-center">
            <div className="row py-3 align-items-center">
                <div className="col-lg-12">
                    {children}
                </div>
            </div>
        </div>
    )
}