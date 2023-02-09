import styles from './styles.module.css';
import React from 'react';

export default function banner (){
    return(
        <a href="https://razomforukraine.org/" className="d4-link-external" style={{
                    'position': 'fixed',
                    'bottom': '0',
                    'left': '-80px',
                    'height': '84px',
                    'width': '300px',
                    'transform': 'rotate(45deg)',
                    'background': 'linear-gradient(-180deg,#1796FC 50%,#ffd500 0)',
                    zIndex: '99999',
                }}></a>
    );
}