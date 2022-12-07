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

def _chess_aq_formulas(formula):
    """Converts a GEMS formula to CHESS format 

    Args:
        formula (string): chemical formula GEMS format

    Returns:
        string: CHESS formula
    """
    if formula == 'H2O@':
        return 'H2O'
    result = formula.find('@')
    if result!=-1:
        return formula.replace("@", "(aq)")
    
    word_list = re.split(r"-", formula)
    if len(word_list)==2:
        return f'{word_list[0]}[{word_list[1]}-]'
    else:
        word_list = re.split(r"\+", formula)
        if len(word_list)==2:
            return f'{word_list[0]}[{word_list[1]}+]'
        else:
            print(f'error in parsing {formula}')

def reaction_to_chess(substance, reaction, Ts, logKs, datasource):
    """Generates a string in CHESS format for a substance with given reaction and data

    Args:
        substance (thermofun::substance): substance object containing data for substance
        reaction (dic): reaction, reactants and their coefficients
        Ts (list): list of temperatures
        logKs (list): list of logKs, should be [0,25,60,100,150,200,250,300]
        datasource (string): source of data

    Returns:
        string: reaction in chess format to be written to txt file
    """
    
    if len(logKs)!=8:
        raise Exception(f"For Chess export the logK list needs to have a size of 8")
        
    TChess = [0,25,60,100,150,200,250,300]
    Ts.sort()
    if Ts!=TChess:
        raise Exception(f"Temperature list is not according to Chess format {TChess}")
    
    try:
        string_composition = f""
        last_key = list(reaction)[-2]
        for x, y in reaction.items():
            if x!=substance.symbol():
                if x == last_key:
                    string_composition = string_composition + f" {y} {to_chess_aq_formulas(x)}"
                else:
                    string_composition = string_composition + f" {y} {to_chess_aq_formulas(x)},"
    
        return (f"  {substance.symbol()} {{\n"
                         f"    vol.weight = {100*substance.molarMass()/(substance.thermoReferenceProperties().volume.val):.2f} kg/m3\n"
                         f"    composition ={string_composition}\n"
                         f"    logK = {logKs[0]:.3f}(0), {logKs[1]:.3f}(25), {logKs[2]:.3f}(60), {logKs[3]:.3f}(100), \\\n"
                         f"           {logKs[4]:.3f}(150), {logKs[5]:.3f}(200), {logKs[6]:.3f}(250), {logKs[7]:.3f}(300)\n"
                         f"    comment {{\n"
                         f"      Ref vol.weight : {datasource}, valid to 100 C\n"
                         f"    }}\n"
                         f"  }}\n")

    except:
        print("An error occurred") 