import React from 'react';
import clsx from 'clsx';
import styles from './styles.module.css';

type FeatureItem = {
  title: string;
  icon: string;
  description: JSX.Element;
  link: string;
};

const FeatureList: FeatureItem[] = [
  {
    title: 'Features',
    icon: 'fa-magic',
    description: (
      <>
        Access, govern, transform, explore, learn, share, apply, deploy and integrate.
      </>
    ),
    link: 'https://datagrok.ai/'
  },
  {
    title: 'Solutions',
    icon: 'fa-briefcase',
    description: (
      <>
        Enterprice, data science, cheminformatics and macromolecules.
      </>
    ),
    link: 'https://datagrok.ai/enterprise'
  },
  {
    title: 'Use cases',
    icon: 'fa-mug-hot',
    description: (
      <>
        Examples of how Datagrok is used by our customers. 
      </>
    ),
    link: 'https://datagrok.ai/solutions'
  },
  {
    title: 'Team',
    icon: 'fa-user-friends',
    description: (
      <>
        Meet the Datagrok team, and see our open positions.
      </>
    ),
    link: 'https://datagrok.ai/company/team'
  },
  {
    title: 'Community',
    icon: 'fa-users-class',
    description: (
      <>
        Ask generic questions, or related to particular entities.
      </>
    ),
    link: 'https://community.datagrok.ai/'
  },
  {
    title: 'Help',
    icon: 'fa-question',
    description: (
      <>
       Explore the help documentation to know how Datagrok works.
      </>
    ),
    link: 'https://datagrok.ai/help/datagrok/'
  },
];

function Feature({title, icon, description, link}: FeatureItem) {
  return (
    <div className='col-6 col-lg-4 col-md-4 col-sm-6 py-3'>
      <div className={"card border-0 shadow-sm p-4 " + styles.featuresCard}>
        <a href={link}>
          <i className={'far fa-lg my-4 '+icon}></i>
          <h3 className='mt-2'>{title}</h3>
          <p style={{color:'var(--grey-5)'}}>{description}</p>
        </a>
      </div>

    </div>
  );
}

export default function HomepageFeatures(): JSX.Element {
  return (
    <section className='bg-light pb-4'>
      <div className="container pb-4 mb-4">
        <div className="row pb-4">
          {FeatureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}
