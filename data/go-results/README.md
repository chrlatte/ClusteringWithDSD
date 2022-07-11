** Format of the GO results
---

The GO results for the DREAM3 network is provided in the `dream_3_go_results.tsv` file.

The file is a tab-separated tsv format. Each field represents
1. GO-LABEL: Name of the Gene Ontology label
2. PROTEIN: The protein name
3. GO-LABEL-TYPE: The type of the protein labels. There are three protein GO labels; `Molecular Function, Biological Process and Cellular Component`

---

Few things to consider:

1. One protein can have many entries in the table (same for GO labels). So this table represents a many to many interaction. If you want to get all the proteins associated with a given GO label, use the following code:

```
import pandas as pd

GO = <Your GO label here>
df = pd.read_csv("dream_3_go_results.tsv", sep = "\t")
proteins = df[df["GO-LABEL"] == GO]["PROTEIN"].values 
### This `proteins` is the required array
```

2. To find all the GO labels associated with a particular protein, use this code block:

```
import pandas as pd

PROT = <Your PROT here>
df = pd.read_csv("dream_3_go_results.tsv", sep = "\t")
gos = df[df["PROTEIN"] == PROT]["GO-LABEL"].values 
### This `gos` is the required array
```

3. For the DREAM3, all the associated GO labels are of level 5 or higher in the GO hierarchy (for all molecular function, biological process and cellular component hierarchies). I only selected the GO labels that labeled at least 10 proteins in the DREAM-3 network.


