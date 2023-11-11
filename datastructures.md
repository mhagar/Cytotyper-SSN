# parsed_dataset (dict)
  ## 1: <cluster>
  ### data (pandas.DataFrame)
  ### gene_key
  ### accession
  ### organism_id
  ### summary (list)
  #### gap
   ##### length (int)
  #### gene_hit
   ##### pfamkey
   ##### direction
   ##### rel_start
   ##### rel_stop
   ##### length
  #### gap

  #### ..
  #### gap
  ### hash (list)
  ## 2: <cluster>
  ## 3: <cluster>
  ## ..
  ## n: <cluster>

# genotypes (dict)
## 1: <genotype>
  ### count int
  ### gene_keys [list]
  ### accessions [list]
  ### organisms [list]
  ### string_hash "str"
## 2: <genotype>
## 3: <genotype>
## ..: <genotype>
## p: <genotype>

# <query>
## data
  - <Dataframe> (Neighbor Table)
## pfams (dict)
  - 1 : [PF4903, PF29032 PF9043..]
  - 2 : [PF42303, PF2352132..]
  - 3 ..
## labels (dict)
  - 1 : 'Prma'
  - 2 : 'Reiske'
  - 3 : 'ATP-grasp'
  - 4 : 'Amidase'

## colors (dict)

## lengths (dict)

# <analysis>
- query
- dataset
- genotypes
- cytotable
