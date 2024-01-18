import Markdown from 'markdown-to-jsx';
import React from 'react';

export default function Cpa({children}){
    return(
        <section className='container-fluid py-5 px-0' style={{backgroundColor: 'var(--brand-accent)'}}>
                <div className="container cpa-container animate__animated animate__fadeIn" >
                    <div class="row  justify-content-center text-center mt-3 mb-5">
                        <div class="col-lg-10">
                            {children}
                        </div>
                </div>
            </div>
        </section>
    )
}