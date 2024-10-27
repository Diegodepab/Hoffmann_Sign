

import pandas as pd
import requests
import stringdb
# cargamos los genes en un DataFrame
genes = pd.read_csv("genes.tsv", sep="\t", header=None, names=["id", "name"])

string_ids = stringdb.get_string_ids(genes["name"], species=9606)
# guardamos las anotaciones en un archivo tsv
string_ids.to_csv("string_ids.tsv", sep="\t", index=False)
print("Data successfully downloaded to string_ids.tsv")
