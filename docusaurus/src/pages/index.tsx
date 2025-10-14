import React from 'react';
import Layout from '@theme/Layout';
import Head from "@docusaurus/Head";
import NavBar from '../components/new/navbar';
import Footer from '../components/new/footer';
import HomeLayout from '../docs/home.mdx';

export default function Home(): JSX.Element {
    return(

                <div className='new'>
                    <NavBar page={''}/>
                    <div class="container-fluid text-center px-0">
                        <HomeLayout/>
                    </div>
                    <Footer/>
                </div>
    );
}
