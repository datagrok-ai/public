import Markdown from 'markdown-to-jsx';
import React from 'react';

const footerIcons = [
    {
        icon: "linkedin-in",
        href: "https://www.linkedin.com/company/datagrok/"
    },
    {
        icon: "youtube",
        href: "https://www.youtube.com/@Datagrok"
    },
    {
        icon: 'github',
        href: 'https://github.com/datagrok-ai'
    }
];

const footerLinks = [
    {
        title: "Solutions",
        links: [
            { name: "Cheminformatics", href: "/cheminformatics" },
            { name: "Bioinformatics", href: "/bioinformatics" },
            { name: "Data Science", href: "/data-science" },
            { name: "Enterprise IT", href: "/enterprise" },
        ],
    },
    {
        title: "Company",
        links: [
            { name: "Careers", href: "/company/careers" },
            { name: "Team", href: "/company/team" },
            { name: "Customer Stories", href: "/solutions" },
        ],
    },
    {
        title: "Resources",
        links: [
            { name: "Wiki", href: "/help" },
            { name: "JS-API", href: "/js-api" },
            { name: "Community", href: "https://community.datagrok.ai" },
            { name: "YouTube", href: "https://www.youtube.com/channel/UCXPHEjOd4gyZ6m6Ji-iOBYg" },
            { name: "User Meetings", href: "https://us02web.zoom.us/meeting/register/up0vfuCrpjgpHNzi371YEJIQ4GkMpTm4disW#/registration" },
        ],
    },
];

function FooterIconsSection({links }) {
    return (
        <div className="ml-auto footer-icons">
            <ul className='p-0'>
                {links.map((link, index) => (
                    <li key={index}>
                        <a href={link.href} target='_blank'><i className={"fa fa-brands fa-" +link.icon}></i></a>
                    </li>
                ))}
            </ul>
        </div>
    )
}

function FooterLinkSection({ title, links }) {
    return (
        <div className="footer-section ml-auto">
            <h6 className='pb-2'>{title}</h6>
            <ul>
                {links.map((link, index) => (
                    <li key={index}>
                        <a href={link.href}>{link.name}</a>
                    </li>
                ))}
            </ul>
        </div>
    )
}

export default function Footer() {
    return (
        <footer className="container-fluid pt-3" style={{backgroundColor: '#1a3248'}}>
            <div className="container pt-5 pb-3">
                <div className="row">
                    <div className='col-lg-4 footer-brand'>
                    <a className="navbar-brand ml-0 mb-0 h1 animate__animated" href="https://datagrok.ai/" > datagrok</a>
                    <p>Our apps are free to use under the MIT license. Everyone is encouraged to add new features and contribute.</p>
                    <a href="mailto:info@datagrok.ai" className='body-link'>info@datagrok.ai</a>
                    </div>
                    <div className='col-lg-8 py-3 d-flex justify-content-end'>
                    {footerLinks.map((section, index) => (
                        <FooterLinkSection key={index} title={section.title} links={section.links} />
                    ))}
                    </div>
                </div>
                <div className="row footer-bottom pt-5"> 
                    <p style={{opacity: "35%"}}>© {new Date().getFullYear()} Datagrok, Inc. All rights reserved.</p>
                    <FooterIconsSection links={footerIcons} /> 
                </div>
            </div>
        </footer>
    )
}
