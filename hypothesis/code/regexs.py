import re
from pathlib import Path

###############################################################################
# Extracting information using Regexps
# 1. Use regular expressions to extract useful information from the cards.
# 2. This information should then be put into the final summary table.

###############################################################################

# Regular expressions to get annotations from a genbank file.
ANNOTATION_REGEXS: list[re.Pattern] = [
    re.compile(r"ACCESSION(.*\d)"),
    re.compile(r"/organism=(\".*\")"),
    re.compile(r"/chromosome=(\".*\")"),
    re.compile(r"/product=(\".*\")"),
    re.compile(r"/gene=(\".*\")"),
    re.compile(r"CDS(.*[0-9])"),
    re.compile(r"DEFINITION(.*|.*\n.*)\.\n"),
]


# ------------------------------------------------------------------------------
# Limpiar str removiendo saltos de linea y comillas.
# ------------------------------------------------------------------------------
def _replace_str(text):
    """
    Clean string removing line breaks and quotes
    param: the text to be cleared
    return: a clean text
    """

    chars = '\n"'
    for c in chars:
        text = text.replace(c, "")
    return text


# Get information with regex.
# ------------------------------------------------------------------------------
def get_info_with_regex(file_path, list_regex: list[re.Pattern]) -> list[str]:
    """
    Get information from genbank with regular expression
    param:  file_path -> genbank file path
            list_regex -> list of required regular expressions
    return: List with the information we request from the genbank
    """
    txt: str = Path(file_path).read_text()
    list_info: list[str] = []
    for pattern in list_regex:

        match_list: list[re.Match] = list(pattern.finditer(txt))

        for match in match_list:
            list_info.append(_replace_str(match.group(1)).strip())

    return list_info
