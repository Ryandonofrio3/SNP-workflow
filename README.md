This program takes in a PLINK file and automatically processes the data, blasts it against AAEG5.0, locates tests and designs primers and finally writes the output to a CSV.

First you will need 2 files:
keeplist.txt & phenotypes.txt

Keeplist.txt
is a two column text file with the samples you are interested in --refer to the Bed file for their exact names:

Example:

| Name     | Samp         |
|----------|--------------|
| LIVPIB12 | LIVPIB12_003 |
| LIVPIB12 | LIVPIB12_002 |
| LIVPIB12 | LIVPIB12_006 |
| LIVPIB12 | LIVPIB12_007 |
| LIVPIB12 | LIVPIB12_008 |
| LIVPIB12 | LIVPIB12_013 |
| LIVPIB12 | LIVPIB12_016 |
| Hanoi    | Han_049      |
| Hanoi    | Han_050      |
| Hanoi    | Han_051      |
| Hanoi    | Han_052      |
| Hanoi    | Han_053      |
| Hanoi    | Han_054      |
| Hanoi    | Han_055      |
| Hanoi    | Han_056      |
| Hanoi    | Han_057      |
| Hanoi    | Han_058      |
| Hanoi    | Han_059      |

Then phenotypes.txt, this is the same as keep list but with a third column
| FID      | IID          | PHENOTYPE |
|----------|--------------|-----------|
| LIVPIB12 | LIVPIB12_003 | 1         |
| LIVPIB12 | LIVPIB12_002 | 1         |
| LIVPIB12 | LIVPIB12_006 | 1         |
| LIVPIB12 | LIVPIB12_008 | 1         |
| LIVPIB12 | LIVPIB12_013 | 1         |
| LIVPIB12 | LIVPIB12_016 | 1         |
| Hanoi    | Han_049      | 2         |
| Hanoi    | Han_050      | 2         |
| Hanoi    | Han_051      | 2         |
| Hanoi    | Han_052      | 2         |
| Hanoi    | Han_053      | 2         |
| Hanoi    | Han_054      | 2         |
| Hanoi    | Han_055      | 2         |
| Hanoi    | Han_056      | 2         |
| Hanoi    | Han_057      | 2         |
| Hanoi    | Han_058      | 2         |
| Hanoi    | Han_059      | 2         |

All you need is a PLINK assoc file. to generate one on the lab comp you can open command prompt and write the following

<br>
<br>

`
cd Desktop/plink2
`

<br>
<br>

`
plink --file file1 --keep keeplist.txt --make-bed --out subset
`
<br>
<br>

`
plink --bfile subset --recode --out recoded_subset
`
<br>
<br>
`
plink --file recoded_subset --pheno phenotypes.txt --assoc --out association_results
`
<br>
<br>
Bring this file into the inputs folder of the Python code and run main.py
The output will be stored in Your name_Date_output.csv


