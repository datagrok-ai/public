import React from 'react';
import clsx from 'clsx';
import styles from './styles.module.css';

export default function heroSection(){
  return (
    <section className="bg-light py-4">
      <div className="container text-center">
        <div className="row">
          <div className="col-12 mt-4">
            <h1 className='mt-4'>Datagrok: Swiss Army Knife for Data</h1>
            <h4 className='mt-2'>A platform for turning data into actionable insights</h4>
          </div>  
        </div>
      </div>
    </section>
  );
}
