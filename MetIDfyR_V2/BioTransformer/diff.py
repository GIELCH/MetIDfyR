import re

def formula_to_dict(formula):
    """Convertit une formule brute en dict : 'C6H4' -> {'C':6, 'H':4}"""
    if not formula:
        return {}
    elements = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    parsed = {}
    for (el, n) in elements:
        parsed[el] = parsed.get(el, 0) + (int(n) if n else 1)
    return parsed

def atom_dict_to_string(atom_dict, priority=['C', 'H', 'N', 'S', 'O', 'F', 'Cl', 'Br', 'P']):
    """
    Transforme un dict d'atomes en une formule chimique unique ordonnée 
    (ex: {'H': 2, 'O': 1} → 'H2O').
    """
    if not atom_dict:
        return ""
    
    # Triez les éléments en respectant la priorité CHNOFS...
    sorted_elements = sorted(atom_dict.items(), 
                             key=lambda x: (priority.index(x[0]) if x[0] in priority else 99, x[0]))
    
    parts = []
    for el, v in sorted_elements:
        if v > 0: # S'assurer que les valeurs non positives sont ignorées
            parts.append(f"{el}{v if v > 1 else ''}")
    return "".join(parts)


def format_formula_diff(diff_dict):
    """
    Transforme un dict de différence en une string formatée où les ajouts 
    sont regroupés (+FORMULE) et les pertes sont regroupées (-FORMULE).
    """
    if not diff_dict:
        return ""

    added = {}
    removed = {}
    
    # 1. Séparation des atomes ajoutés et retirés
    for el, delta in diff_dict.items():
        if delta > 0:
            added[el] = delta
        elif delta < 0:
            removed[el] = abs(delta)

    # 2. Conversion des dicts en formules chimiques uniques
    formula_removed = atom_dict_to_string(removed)
    formula_added = atom_dict_to_string(added)

    output = []

    # 3. Formatage de la sortie (-FORMULE +FORMULE)
    if formula_removed:
        output.append(f"-{formula_removed}")

    if formula_added:
        output.append(f"+{formula_added}")

    # 4. Nettoyage final pour le style
    return " ".join(output).strip()


def diff_formulas(form1, form2):
    """Calcule la différence entre deux formules brutes (form2 - form1)."""
    dict1 = formula_to_dict(form1)
    dict2 = formula_to_dict(form2)
    
    all_elements = set(dict1.keys()) | set(dict2.keys())
    
    diff = {}
    for el in all_elements:
        delta = dict2.get(el, 0) - dict1.get(el, 0)
        if delta != 0:
            diff[el] = delta
            
    return format_formula_diff(diff)
  
def smiles_to_image(smiles, outfile):
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Draw.MolToFile(mol, outfile, size=(400, 400))
    return outfile
