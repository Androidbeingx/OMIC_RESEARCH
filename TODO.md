# TODO

## Choose the following:
- An organism A.
- Several species of the organism A.
- At least one G gene present in the species.

Download G gene information for each species, extract gene information, align the sequences, and try to write some conclusions.
- The important thing is not the conclusion, but the process carried out.
- Approximately 10 sequences must be compared.
- If you do not have enough species to reach all 10 sequences, choose more than one organism.
- Each group must study a different gene and organism.

1. Formulate the hypothesis (1.00 points)
- Formulate a hypothesis that you want to validate about the gene and the chosen species.
- The important thing is to do it in a logical and meaningful way.
- It does not matter its scientific relevance or whether the answer is already known or not.

2. Download Genbanks using Bio.Entrez (1.00 points)
- Search for information using Bio.Entrez and choose the correct fields. You can pre-test manually in the NCBI web interface.
- Use BioPython to automatically download GenBank files.

3. Extract information using Regexps (1.00 points)
- Genbank files contain a lot of textual information.
- Use regular expressions to extract useful information from records.
- You will then have to put this information in the final summary table.

4. Sequence alignment (1.00 points)
- Use BioPython to align one of the sequences with the other sequences.
- There must be at least 10 comparisons.
- Save the alignment score and order the comparisons according to similarity.

5. Show the results (1.00 points)
- Show a summary table made in Pandas showing at least the following columns:
  - Accession number of the sequence.
  - Descriptive name of the sequence.
  - Relevant information extracted with regexps.
  - Alignment score.
- The table must be ordered by the score of the alignment, from highest to lowest score.
- Study the similarity of the highest scoring sequences and their Genbank tokens.
- Try to draw some conclusion and validate or not your hypothesis. A paragraph or two is enough.
- The important thing is not the conclusion but to see that you have done the whole process correctly.
- The code, table, and conclusions must be in a JupyterLab notebook.

## Evaluation
- All students must have programmed Python code.
- Code must work to be evaluated, otherwise the practice will be scored with a zero.
- The students will have to do a demonstration and answer the teacher's questions. The teacher can request changes to the code.
- Students must correctly answer the questions about their own code and make the requested changes (if any), otherwise the practice will be scored with a zero.

## Resources
- NCBI: https://www.ncbi.nlm.nih.gov/
- BioPython: https://biopython.org/
- EOL: https://eol.org/
