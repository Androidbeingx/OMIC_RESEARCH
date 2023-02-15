import re
from pathlib import Path

###############################################################################
    # Extracting information using Regexps
    # 1. Use regular expressions to extract useful information from the cards.
    # 2. This information should then be put into the final summary table.

###############################################################################

# Path file genbak
FILE_PATH = Path('output/alcohol_dehydrogenase-drosophila_melanogaster.gb')

# Regexps
gene:       str     = r'gene(.*[0-9]\n)'
cds:        str     = r'CDS(.*[0-9])'
product:    str     = r'/product=(\".*\")'
mol_type:   str     = r'/mol_type=(\".*\")'
chromosome: str     = r'/chromosome=(\".*\")'
organism:   str     = r'/organism=(\".*\")'
definition: str     = r'DEFINITION(.*\n.*)'
accession:  str     = r'ACCESSION(.*\d)'

#------------------------------------------------------------------------------
def get_info_with_regex(filename, regex,key):

    txt:    str     =  Path(filename).read_text()
    reg:    str     = regex
    pat:re.Pattern  = re.compile(reg)
    
    match_list: list[re.Match] = list(pat.finditer(txt))

    dict_info: dict = {}
    for match in match_list:
        val = match.group(1).replace('\n','').strip()
        dict_info = {key: val}

    print(dict_info)
    #return dict_info

#------------------------------------------------------------------------------
def main():
    #read_genbank_file('output/alcohol_dehydrogenase-drosophila_melanogaster.gb')
    get_info_with_regex(FILE_PATH, gene, "gene")
    get_info_with_regex(FILE_PATH, cds, "cds")
    get_info_with_regex(FILE_PATH, product,"product")
    get_info_with_regex(FILE_PATH, mol_type, "mol_type")
    get_info_with_regex(FILE_PATH, chromosome, "chromosome")
    get_info_with_regex(FILE_PATH, organism, "organism")
    get_info_with_regex(FILE_PATH, definition, "definition")
    get_info_with_regex(FILE_PATH, accession, "accession")

#------------------------------------------------------------------------------
# Main
if __name__ == "__main__":
    main()
