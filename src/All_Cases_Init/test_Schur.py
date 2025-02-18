import subprocess


import subprocess
import re

def run_schur(commands):
    """Exécute Schur et filtre la sortie pour ne conserver que le résultat du pléthysme."""
    process = subprocess.run(['schur'], input=commands, text=True, capture_output=True)
    output_lines = process.stdout.split("\n")  # Découper en lignes
    
    # Filtrage des lignes contenant le résultat des partitions (on cherche "{")
    result = ""
    capture = False
    for line in output_lines:
        #if "pl " in line:  # Début de la partie utile
        #    capture = True
        if "{" in line:  # Lignes contenant les partitions
            result+=line.strip()

    ## Reading partitions and coefficients
    partition_dict = {}

    # Regex pour capturer les termes {partition} ou coef{partition}
    matches = re.findall(r'(\d*)\{([\d\s\^]+)\}', result)
    for coef, partition in matches:
        coef=int(coef) if coef else 1
        partition_expanded = re.sub(r'(\d+)\^(\d+)', lambda m: (m.group(1) + " ") * int(m.group(2)), partition).strip() # Repeat k times when ^k
        partition_list = tuple(map(int, partition_expanded.split()))
        partition_dict[partition_list]=coef
        # TODO : gérer les nombres à deux chiffres avec ! dans Schur. Voir le manuel ici :
        # https://master.dl.sourceforge.net/project/schur/manual/6.07/schur_manual.pdf?viasf=1
        # Exemple : !15!13!11 975 signifie [15,13,11,9,7,5]
        # Exemple : !137!29!10 941 signifie [137,29,10,9,4,1] 
        # Exemple : !12ˆ!13 !10ˆ!31 9ˆ3 21 would correspond to the partition [12]*13+[10]*31+[9]*3+[2,1].
    print(result)    
    return partition_dict  # Reformater proprement

# Commandes à exécuter
commands = """gr u4
rep
pl 321,21
"""

# Exécution et affichage du résultat propre
output = run_schur(commands)
print(output)
