import re

def conversion(path):
    # Lire le fichier .m
    with open(path, 'r') as file:
        content = file.read()

    # Fonction pour extraire les données sous forme de dictionnaire
    def extract_matrix_to_dict(text, headers, costtype=False):
        matrix_pattern = re.compile(r'\[([^\]]+)\]', re.DOTALL)
        match = matrix_pattern.search(text)
        if match:
            matrix_text = match.group(1).strip()
            rows = matrix_text.split(';')
            data = []
            for row in rows:
                values = list(map(float, row.split()))
                if costtype:
                    if(len(values)>1):
                        # For gencost, format it in the specific dictionary
                        cost_dict = {
                            'costtype': values[0],  # First value is costtype
                            'startup': values[1],   # Second value is startup
                            'shutdown': values[2],  # Third value is shutdown
                            'n': values[3],         # Fourth value is n
                            'costvector': values[4:]  # Remaining values form costvector
                        }
                        data.append(cost_dict)
                else:
                    data.append(dict(zip(headers, values)))
            return data
        return None

    # Dictionnaire pour stocker toutes les données
    mpc = {}

    # Extraire le nom du cas (casename)
    casename_match = re.match(r'function\s+mpc\s*=\s*(\w+)', content)
    if casename_match:
        mpc['casename'] = casename_match.group(1)

    # Extraire la donnée des buses (bus data)
    bus_data_text = re.search(r'mpc\.bus\s*=\s*\[.*?\];', content, re.DOTALL)
    if bus_data_text:
        bus_headers = [
            'bus_i', 'type', 'Pd', 'Qd', 'Gs', 'Bs', 'area', 'Vm', 'Va', 'baseKV', 'zone', 'Vmax', 'Vmin'
        ]
        mpc['bus'] = extract_matrix_to_dict(bus_data_text.group(0), bus_headers)
    if(mpc['bus'][-1]=={}):
        mpc['bus']=mpc['bus'][:-1]


    # Extraire les données des générateurs (generator data)
    gen_data_text = re.search(r'mpc\.gen\s*=\s*\[.*?\];', content, re.DOTALL)
    if gen_data_text:
        gen_headers = [
            'bus', 'Pg', 'Qg', 'Qmax', 'Qmin', 'Vg', 'mBase', 'status', 'Pmax', 'Pmin', 'Pc1', 'Pc2', 'Qc1min', 
            'Qc1max', 'Qc2min', 'Qc2max', 'ramp_agc', 'ramp_10', 'ramp_30', 'ramp_q', 'apf'
        ]
        mpc['gen'] = extract_matrix_to_dict(gen_data_text.group(0), gen_headers)
    if(mpc['gen'][-1]=={}):
        mpc['gen']=mpc['gen'][:-1]

    # Extraire les données des branches (branch data)
    branch_data_text = re.search(r'mpc\.branch\s*=\s*\[.*?\];', content, re.DOTALL)
    if branch_data_text:
        branch_headers = [
            'fbus', 'tbus', 'r', 'x', 'b', 'rateA', 'rateB', 'rateC', 'ratio', 'angle', 'status', 
            'angmin', 'angmax'
        ]
        mpc['branch'] = extract_matrix_to_dict(branch_data_text.group(0), branch_headers)
    if(mpc['branch'][-1]=={}):
        mpc['branch']=mpc['branch'][:-1]
    for branch in range(len(mpc['branch'])):
        mpc['branch'][branch]["fbus"]=int(mpc['branch'][branch]["fbus"])
        mpc['branch'][branch]["tbus"]=int(mpc['branch'][branch]["tbus"])

    # Extraire la donnée de base MVA
    baseMVA_text = re.search(r'mpc\.baseMVA\s*=\s*(\d+);', content)
    if baseMVA_text:
        mpc['baseMVA'] = float(baseMVA_text.group(1))

    # Extraire la donnée de cost (gencost)
    gencost_text = re.search(r'mpc\.gencost\s*=\s*\[.*?\];', content, re.DOTALL)
    if gencost_text:
        # Extract gencost and format it
        mpc['gencost'] = extract_matrix_to_dict(gencost_text.group(0), None, costtype=True)
    if(mpc['gencost'][-1]=={}):
        mpc['gencost']=mpc['gencost'][:-1]

    # Afficher les données extraites
    return mpc
