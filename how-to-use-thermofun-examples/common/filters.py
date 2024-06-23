# Copyright (C) 2022 dmiron
# 
# This file is part of thermoexport.
# 
# thermoexport is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# thermoexport is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with thermoexport.  If not, see <http://www.gnu.org/licenses/>.

import json
import jmespath as jp

def substances_containing_elements(dataset, matches, all_=True):
    """prints symbols of substances containing elements 
    Args:
        dataset (string): thermofun dataset in json string
        matches (list): elements to match
        all_ (bool, optional): match all or any from the elements list. Defaults to True.
    """

    # Opening JSON file
    f = open(dataset)
    data = json.load(f)

    results = jp.search("substances[*].{formula: formula, symbol: symbol, name: name, class_: class_.*}", data)

    print("{:<20} {:<40} {:<20}".format('symbol','formula', 'class'))

    for r in results:
        if all_:
            if all(x in r['formula'] for x in matches):
                print("{:<20} {:<40} {:<20}".format(r['symbol'],r['formula'], r['class_'][0]))
                #print(f'symbol \'{r['symbol']}\' ')
        else:
            if any(x in r['formula'] for x in matches):
                print("{:<20} {:<40} {:<20}".format(r['symbol'],r['formula'], r['class_'][0]))
                #print(r['symbol'])
                

def solids_containing_elements(dataset, matches, all_=True):
    """prints symbols of substances containing elements 
    Args:
        dataset (string): thermofun dataset in json string
        matches (list): elements to match
        all_ (bool, optional): match all or any from the elements list. Defaults to True.
    """

    # Opening JSON file
    f = open(dataset)
    data = json.load(f)

    results = jp.search("substances[*].{formula: formula, symbol: symbol, name: name, class_: class_.*}", data)

    print("{:<20} {:<40} {:<20}".format('symbol','formula', 'class'))

    for r in results:
        if all_:
            if all(x in r['formula'] for x in matches):
                if r['class_'][0] == '':
                    print("{:<20} {:<40} {:<20}".format(r['symbol'],r['formula'], r['class_'][0]))
                #print(f'symbol \'{r['symbol']}\' ')
        else:
            if any(x in r['formula'] for x in matches):
                if r['class_'][0] == '':
                    print("{:<20} {:<40} {:<20}".format(r['symbol'],r['formula'], r['class_'][0]))
                #print(r['symbol'])