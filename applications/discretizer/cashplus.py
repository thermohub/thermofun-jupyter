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

import xgems as xg
import numpy as np
import pandas as pd
import itertools
import json
import chemicalfun as cf
from common.pyGEMS import GEMS, GEMSParameter

def from_ratio_to_formula(ratio):
    if (ratio=='Ca/Si'):
        return {'Ca':1,'O':1}
    elif (ratio=='Na/Si'):
        return {'Na':1,'O':0.5}
    elif (ratio=='K/Si'):
        return {'K':1,'O':0.5}
    elif (ratio=='Al/Si'):
        return {'Al':1,'O':1.5}

def make_list_of_lists(setup_):
    listOLists = [ c['values'] for c in setup_]
    results = list(itertools.product(*listOLists))
    return results


def discretize(composition):
    setup_ = [{'formula': from_ratio_to_formula(c['ratio']), 'values': [ x for x in np.arange(c['interval'][0], c['interval'][1]+1e-11, c['step']) ]} for c in composition]
    
    discretized_pairs_list = make_list_of_lists(setup_)
    
    result = [[{'formula': setup_[i]['formula'], 'value': v} for i, v in enumerate(pair)] for pair in discretized_pairs_list]
        
    return result

def from_CASH_elements_to_oxides(composition):
    Si= composition['Si']
    Ca= composition['Ca']
    Na= composition['Na']
    K=  composition['K']
    Al= composition['Al']
    H= composition['H']
    
    Ox = 0.0
    result = []
    
    if Si>1e-14:
        result.append({'SiO2': Si})
        Ox = Ox + (2/1)*Si
    if Ca>1e-14:
        result.append({'CaO': Ca})
        Ox = Ox + (1/1)*Ca
    if Na>1e-14:
        result.append({'Na2O': 0.5*Na})
        Ox = Ox + (1/2)*Na
    if K>1e-14:
        result.append({'K2O': 0.5*K})
        Ox = Ox + (1/2)*K
    if Al>1e-14:
        result.append({'Al2O3': 0.5*Al})
        Ox = Ox + (3/2)*Al
    if H>1e-14:
        result.append({'H2O': 0.5*H})
        Ox = Ox + (1/2)*H
    
    if abs((Ox-composition['O']))>1e-14:
        print (Ox-composition['O'])

    return result

def calculate_discretized_model(composition, discretized_composition, digits=4):
    #path to the GEMS system file, change this if a different gems system file is used
    gems_system_filename = 'discretizer/gems/CASH+discr-dat.lst' 

    gems = GEMS(gems_system_filename)        # initialize GEMS with the system file
    
    gems.clear()
    
    cash = []
    
    for d in discretized_composition:
        gems.clear()
        gems.add_amt_from_formula({'O':2},         0.001, 'moles')
        gems.add_amt_from_formula({'Nit':0.5},         0.008, 'moles')
        gems.add_amt_from_formula({'H':2,'O':1},         2, 'moles')
        gems.add_amt_from_formula({'Si':1,'O':2},         1, 'moles')
        s = f'For '
        ratio = 'ratio'
        value = 'value'
        for i,c in enumerate(d):
            s = s + f'{composition[i][ratio]} {round(c[value], digits)} '
            gems.add_amt_from_formula(c['formula'], c['value'], 'moles')
        s = s + f'{gems.equilibrate()}\n'
        
        cash_composition = from_CASH_elements_to_oxides(gems.phases_elements_moles['CASH+'])
        cash_mg = gems.phases_molar_gibbs_energy['CASH+'] # J/mol
        cash_mh = gems.phases_molar_enthalpy['CASH+']
        cash_ms = gems.phases_molar_entropy['CASH+']
        cash_mcp = gems.phases_molar_heat_capacity['CASH+']
        cash_mv = gems.phases_molar_volume['CASH+'] # cm3
        cash_amount = gems.phases_moles['CASH+'] # moles
        cash_mm = gems.phases_mass['CASH+']/cash_amount # g/mol
        sio2='SiO2'
        cao='CaO'
        formula = f''
        for c in cash_composition:
            for k in c:
                formula = formula + f'({k}){round(c[k], digits)}'
        
        s = s+formula
        print(s)
        
        # cash mg gibbs energy is divded with amount as workaround to a big in gems3k api
        cash.append({'composition': cash_composition, 'formula': formula, 'G': cash_mg, 'H':cash_mh*cash_amount, 'S': cash_ms*cash_amount,
                    'Cp': cash_mcp*cash_amount, 'V': cash_mv*1e6*cash_amount, 'Mmas': cash_mm*1000*cash_amount})
    return cash

def make_fun_cash_phase(cash_phase, j):
    
    j['symbol']=cash_phase['formula']
    j['name']=cash_phase['formula']
    j['formula']=cash_phase['formula']
    
    j['TPMethods'][0]['m_heat_capacity_ft_coeffs']['values'][0]=cash_phase['Cp']
    
    j['mass_per_mole']=cash_phase['Mmas']
    j['sm_gibbs_energy']['values'][0]=cash_phase['G']
    j['sm_enthalpy']['values'][0]=cash_phase['H']
    j['sm_entropy_abs']['values'][0]=cash_phase['S']
    j['sm_heat_capacity_p']['values'][0]=cash_phase['Cp']
    j['sm_volume']['values'][0]=cash_phase['V']/10
    
    j['datasources'][0] = 'discretizer'
    
    
    return j
    
def add_cash_to_fundb(fun_db, cash):
    
    for c in cash:
        #fun_cemdata18.addSubstance('''
        #{"Pst": 100000, "TPMethods": [{"limitsTP": {"lowerT": 273.15, "range": true, "upperT": 373.15}, "m_heat_capacity_ft_coeffs": {"names": ["a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "a10"], "units": ["J/(mol*K)", "J/(mol*K^2)", "(J*K)/mol", "J/(mol*K^0.5)", "J/(mol*K^3)", "J/(mol*K^4)", "J/(mol*K^5)", "(J*K^2)/mol", "J/mol", "J/(mol*K^1.5)", "J/(mol*K)"], "values": [348.02928015774563, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]}, "method": {"0": "cp_ft_equation"}}, {"method": {"34": "mv_constant"}}], "Tst": 298.15, "aggregate_state": {"3": "AS_CRYSTAL"}, "class_": {"0": "SC_COMPONENT"}, "datasources": ["discretizer"], "formula": "(SiO2)1.0(CaO)0.7(H2O)1.0391", "formula_charge": 0, "mass_per_mole": 348.1544652586432, "name": "(SiO2)1.0(CaO)0.7(H2O)1.0391", "sm_enthalpy": {"units": ["J/mol"], "values": [-5050223.333329548]}, "sm_entropy_abs": {"units": ["J/(mol*K)"], "values": [369.32087843116216]}, "sm_gibbs_energy": {"units": ["J/mol"], "values": [-4664977.116092457]}, "sm_heat_capacity_p": {"units": ["J/(mol*K)"], "values": [348.02928015774563]}, "sm_volume": {"units": ["J/bar"], "values": [13.856891906248373]}, "symbol": "(SiO2)1.0(CaO)0.7(H2O)1.0391"}''')
        #sdata = f'{json.dumps(make_fun_cash_phase(c, template))}'
        #print(sdata)
        fun_cemdata18.addSubstance(json.dumps(make_fun_cash_phase(c, template)))
        #print(c['formula'])
        
        
    return fun_cemdata18

def generate_reactions_one_by_one(masters, products, formation_=True):
    reaction_dic = []
    reaction_list = []
    for p in products:
        chemicalReactions = cf.ChemicalReactions(masters+[p])
        reactions = chemicalReactions.generateReactions(formation=formation_) # returns the reactions list as a list of tuples ('substance', coefficient)

        # list of dictionaries, with reaction substances as keys and the reaction coefficients as values
        reaction_dic = reaction_dic+[{el[0]: el[1] for el in r} for r in reactions] 

        # strings of reactions that can be used in ThermoFun to calculate the logK at different T and P
        reaction_list = reaction_list+chemicalReactions.stringReactions()
    for r in reaction_list:
        print(f'{r}')
    return reaction_dic, reaction_list


def discretize_cash(discretization, digits=4):
    discretized_composition = discretize(discretization)
    
    return calculate_discretized_model(discretization, discretized_composition, digits)

def make_fun_cash_phase(cash_phase, j):
    
    j['symbol']=cash_phase['formula']
    j['name']=cash_phase['formula']
    j['formula']=cash_phase['formula']
    
    j['TPMethods'][0]['m_heat_capacity_ft_coeffs']['values'][0]=cash_phase['Cp']
    
    j['mass_per_mole']=cash_phase['Mmas']
    j['sm_gibbs_energy']['values'][0]=cash_phase['G']
    j['sm_enthalpy']['values'][0]=cash_phase['H']
    j['sm_entropy_abs']['values'][0]=cash_phase['S']
    j['sm_heat_capacity_p']['values'][0]=cash_phase['Cp']
    j['sm_volume']['values'][0]=cash_phase['V']/10
    
    j['datasources'][0] = 'discretizer'
    
    
    return j
    
def add_cash_to_fundb(fun_db, cash):

    # Opening JSON file
    f = open('discretizer/template_cash_solid.json')
  
    # returns JSON object as 
    # a dictionary
    template = json.load(f)
    
    for c in cash:
        #fun_cemdata18.addSubstance('''
        #{"Pst": 100000, "TPMethods": [{"limitsTP": {"lowerT": 273.15, "range": true, "upperT": 373.15}, "m_heat_capacity_ft_coeffs": {"names": ["a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "a10"], "units": ["J/(mol*K)", "J/(mol*K^2)", "(J*K)/mol", "J/(mol*K^0.5)", "J/(mol*K^3)", "J/(mol*K^4)", "J/(mol*K^5)", "(J*K^2)/mol", "J/mol", "J/(mol*K^1.5)", "J/(mol*K)"], "values": [348.02928015774563, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]}, "method": {"0": "cp_ft_equation"}}, {"method": {"34": "mv_constant"}}], "Tst": 298.15, "aggregate_state": {"3": "AS_CRYSTAL"}, "class_": {"0": "SC_COMPONENT"}, "datasources": ["discretizer"], "formula": "(SiO2)1.0(CaO)0.7(H2O)1.0391", "formula_charge": 0, "mass_per_mole": 348.1544652586432, "name": "(SiO2)1.0(CaO)0.7(H2O)1.0391", "sm_enthalpy": {"units": ["J/mol"], "values": [-5050223.333329548]}, "sm_entropy_abs": {"units": ["J/(mol*K)"], "values": [369.32087843116216]}, "sm_gibbs_energy": {"units": ["J/mol"], "values": [-4664977.116092457]}, "sm_heat_capacity_p": {"units": ["J/(mol*K)"], "values": [348.02928015774563]}, "sm_volume": {"units": ["J/bar"], "values": [13.856891906248373]}, "symbol": "(SiO2)1.0(CaO)0.7(H2O)1.0391"}''')
        #sdata = f'{json.dumps(make_fun_cash_phase(c, template))}'
        #print(sdata)
        fun_db.addSubstance(json.dumps(make_fun_cash_phase(c, template)))
        #print(c['formula'])
        
        
    return fun_db