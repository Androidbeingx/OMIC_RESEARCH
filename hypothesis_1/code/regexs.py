import re
from pathlib import Path

###############################################################################
    # Extracting information using Regexps
    # 1. Use regular expressions to extract useful information from the cards.
    # 2. This information should then be put into the final summary table.

###############################################################################

# PATH FILE GENBANK
# list_files: list[Path] = ['dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_ananassae.gb',
#                           'dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_erecta.gb',
#                           'dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_melanogaster.gb']


# REGEX
list_regex: list[str] = [
    r'ACCESSION(.*\d)',
    r'/organism=(\".*\")',
    r'/chromosome=(\".*\")',
    r'/product=(\".*\")',
    r'/gene=(\".*\")',
    r'CDS(.*[0-9])',
    r'DEFINITION(.*|.*\n.*)\.\n'
]


#------------------------------------------------------------------------------
# Limpiar str removiendo saltos de linea y comillas. 
#------------------------------------------------------------------------------
def replace_str(text):
    """
    Clean string removing line breaks and quotes
    param: the text to be cleared
    return: a clean text
    """

    chars = '\n"'
    for c in chars:
        text = text.replace(c,'')
    return text


# Get information with regex.
#------------------------------------------------------------------------------
def get_info_with_regex(file_path, list_regex)-> list[str]:
    """
    Get information from genbank with regular expression
    param:  file_path -> genbank file path
            list_regex -> list of required regular expressions
    return: List with the information we request from the genbank
    """
    txt:    str     =  Path(file_path).read_text()
    list_info:list[str] =[]
    for reg in list_regex:
        pat:re.Pattern  = re.compile(reg)
    
        match_list: list[re.Match] = list(pat.finditer(txt))

        for match in match_list:
            list_info.append(replace_str(match.group(1)).strip())
        
    return list_info

#------------------------------------------------------------------------------
def create_dict_info_gb(file_path)-> dict[str,str]:
    """
    Create dictionary to identify what each value is.
    param: file_path -> the file path.
    return: dict_info_gb -> information organized with your key.
    """

    list_key: list[str] = [
        'accession',
        'organism',
        'chromosome',
        'product',
        'gene',
        'cds',
        'definition']
    
    list_value: list[str] = get_info_with_regex(file_path, list_regex)

    dict_info_gb = dict(zip(list_key,list_value))
    return dict_info_gb

#------------------------------------------------------------------------------
def get_info_all_files_gb(list_files)-> list[dict]:
    """
    Get information from a list of gb files.
    param: list_files -> list of file paths where the information is sought.
    return: list_dict_info -> list of dictionaries with information.
    """

    list_dict_info:list[dict] =[]
    for file_path in list_files:
        list_dict_info.append(create_dict_info_gb(file_path))

    return list_dict_info

#------------------------------------------------------------------------------
def main():

    #info = get_info_all_files_gb(list_files)
    #print(info)
    return ""

#------------------------------------------------------------------------------
# Main
if __name__ == "__main__":
    main()