/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/* eslint-disable global-require */

import React from 'react';
import clsx from 'clsx';

interface Props {
  name: string;
  url: string;
  description: string[];
  button_text: string,
}

export default function Card({name, description, url, button_text}: Props) {
  return (
    <div className="col col--6 margin-bottom--lg" style={{'max-width': '350px'}}>
      <div className={clsx('card')}>
        <div className="card__body">
          <h3 style={{'marginTop': '0px'}}>{name}</h3>
          <ul style={{'paddingLeft': '15px'}}>
          {
            description.map((point) => (<li> {point} </li>))
          }
          </ul>
        </div>
        <div className="card__footer">
          <div className="button-group button-group--block">
            <a href={url} className="button button--secondary">
                {button_text}
            </a>
          </div>
        </div>
      </div>
    </div>
  );
}