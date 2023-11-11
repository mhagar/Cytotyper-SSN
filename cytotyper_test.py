import pandas as pd
import sqlite3
import os
from statistics import mode

gap_unit = 500  # Number of bps before adding a 'gap' in the BGC legend
min_occurence = 4  # Number of times a BGC type has to occur or it'llget tossed
window_width = 10  # Number of 'steps' front and back to consider (center=stzf)

# See "datastructures.md" for a diagram of the datastructures.


def importdata(rawdata):
    # connect to sqlite file
    conn = sqlite3.connect(os.path.join('uploads', rawdata.filename))

    sqlquery = "id, num, family, ipro_family, rel_start, rel_stop, \
                direction, gene_key, start, color"
    main_df = pd.read_sql_query(f"SELECT {sqlquery} FROM neighbors", conn)

    # 'attributes' table in db, index corresponds to gene key in 'neighbors'!
    # i.e. gene_key 573, which refers to kutzneria bgc, the stzF homolog is
    # found in 'attributes' row 573.

    sqlquery = "sort_key, accession, num, family,\
                 ipro_family, direction, start, strain"
    organism_df = pd.read_sql_query(f"SELECT {sqlquery} FROM attributes", conn,
                                    index_col='sort_key')

    # 'AHMAD-BORSCHT-CAR' -> ['AHMAD', 'BORSCHT', 'CAR']
    def split_dash(x):
        x = x.split('-')
        return x

    main_df['family'] = main_df['family'].apply(split_dash)
    main_df['ipro_family'] = main_df['ipro_family'].apply(split_dash)

    # combine 'family' and 'ipro_family' columns
    main_df['family'] += main_df['ipro_family']
    main_df.drop(labels='ipro_family', axis=1, inplace=True)
    return main_df, organism_df


class query():
    # Just a test query for now.
    def __init__(self, pfams, labels):
        self.data = 'call parse_dataset() first'
        self.pfams = pfams
        self.labels = labels
        # For now this generates an empty dict {1: [], 2: []}
        # for each item in query.pfams
        self.colors = {k: [] for k, v in self.pfams.items()}
        self.lengths = {k: [] for k, v in self.pfams.items()}

    def collapse_colors_lengths(self):
        self.colors = {k: mode(v) for k, v in self.colors.items()}
        self.lengths = {k: mode(v) for k, v in self.lengths.items()}

        # in case the mode returns multiple values instead of one:
        self.colors = {k: v.split(',')[0] for k, v in self.colors.items()}
        # self.lengths = {k: v.split(',')[0] for k, v in self.lengths.items()}


class gene_hit:
    def __init__(self, pfamkey, direction,
                 rel_start, rel_stop):
        self.pfamkey = pfamkey
        self.direction = direction
        self.rel_start = rel_start
        self.rel_stop = rel_stop
        self.length = abs(rel_start - rel_stop)


class gap:
    def __init__(self, length):
        self.length = length


def groupdata(df, gene_key):
    # 'grouped_df.get_group(573)' returns a df with only the rows with
    #  genekey = 573
    if gene_key in df['gene_key'].values:
        grouped_df = df.groupby(by="gene_key")
        return grouped_df.get_group(gene_key).reset_index()
    else:
        return 0


class cluster:
    def __init__(self, data, gene_key, accession, organism_id):
        self.data = data
        self.gene_key = gene_key
        self.accession = accession
        self.organism_id = organism_id
        self.summary = []
        self.hash = []


def parse_cluster(cluster, attributes, query):
    # for every row in the gene cluster table..
    data = cluster.data
    if isinstance(data, pd.DataFrame) is False:
        return "Null"

    # Generating cluster.summary:
    for i in range(0, len(data)):
        hit_flag = 0

        # if you just passed the row where "stzF" would be..
        if data.loc[i, 'start'] > attributes['start']:
            # AND if you didnt already place "stzF" in this cluster.summary..
            if 0 not in [x.pfamkey for x in cluster.summary
                         if isinstance(x, gene_hit)]:
                # add an entry for "stzF".
                hit = gene_hit(pfamkey=0,
                               direction=attributes['direction'],
                               rel_start=0,
                               rel_stop=0)
                cluster.summary.append(hit)

        # for every item in the query .
        for pfamkey, pfams in query.pfams.items():
            # if there's a match between a label's pfams and the cluster table-
            if [x for x in pfams if x in data.loc[i, 'family']] != []:
                # create a 'gene_hit' object..
                hit = gene_hit(pfamkey=pfamkey,
                               direction=data.loc[i, 'direction'],
                               rel_start=data.loc[i, 'rel_start'],
                               rel_stop=data.loc[i, 'rel_stop'])

                cluster.summary.append(hit)
                hit_flag = 1

                # append the hit's color + length to query.color and .length
                query.colors[pfamkey].append(data.loc[i, 'color'])
                query.lengths[pfamkey].append(hit.length)

        # if this row didn't match anything in the query..
        if hit_flag == 0:
            gap_length = abs(data.loc[i, 'rel_start']
                             - data.loc[i, 'rel_stop'])
            # if there's already a gap object at the end of clusterlist..
            if len(cluster.summary) > 0 \
                    and isinstance(cluster.summary[-1], gap):
                # extend its length property.
                cluster.summary[-1].length += gap_length
            else:
                # add a gap object.
                cluster.summary.append(gap(gap_length))

    # Generating cluster.hash:
    for k, item in enumerate(cluster.summary):
        if isinstance(item, gap):
            # If you're at the beginning OR end of the cluster.summary..
            if k == 0 or k+1 == len(cluster.summary):
                # Append ONE 'G'
                cluster.hash.append('G')
            else:
                # Otherwise, append a 'G' for every 'gap_unit' basepairs
                # default: 300
                for n in range(0, (item.length // gap_unit)):
                    cluster.hash.append('G')

        if isinstance(item, gene_hit):
            direction = 1
            if item.direction == 'complement':
                direction = -1
            cluster.hash.append(str(item.pfamkey * direction))

    # Trimming cluster.hash to the window_width variable:
    if '0' in cluster.hash:
        center = cluster.hash.index('0')
        start = max((center - window_width), 0)  # subtract but minimum 0
        end = min((center + window_width), (len(cluster.hash) - 1))
        cluster.hash = cluster.hash[start:end]

        # removing all trailing 'G's at beginning and end
        m = {k: v for (k, v) in enumerate(cluster.hash) if v != 'G'}
        start, end = min(m.keys()), (max(m.keys()) + 1)
        cluster.hash = cluster.hash[start:end]


def parse_dataset(neighbor_table, attributes_table, query):
    parsed_dataset = {}
    data = attributes_table
    query.data = neighbor_table
    for i in range(1, len(data)+1):
        c = cluster(data=groupdata(neighbor_table, i),
                    gene_key=i,
                    accession=data['accession'].loc[i],
                    organism_id=data['strain'].loc[i])
        parse_cluster(c, data.loc[i], query)
        if '0' in c.hash:  # for some reason some clusters dont have stzF (??)
            parsed_dataset[i] = c

    # collapse the list of colors+lengths encountered into single values
    query.collapse_colors_lengths()

    return parsed_dataset


class genotype:
    def __init__(self, string_hash):
        self.count = 0
        self.gene_keys = []
        self.accessions = []
        self.organisms = []
        self.string_hash = string_hash


def bin_hashes(dataset):

    # Go through all the clusters and create a new dict entry for every new
    # hash encountered. If it's already in the dict, increment the count by one
    unique_string_hashes = {}
    genotypes = {}

    for key, cluster in dataset.items():
        string_hash = ",".join(cluster.hash)

        # if this string-hash hasn't been encountered before..
        if string_hash not in unique_string_hashes:
            # create a new entry
            unique_string_hashes[string_hash] = 0
            genotypes[string_hash] = genotype(string_hash=string_hash.split(
                                              ","))

        unique_string_hashes[string_hash] += 1
        genotypes[string_hash].count += 1
        genotypes[string_hash].gene_keys.append(cluster.gene_key)
        genotypes[string_hash].accessions.append(cluster.accession)
        genotypes[string_hash].organisms.append(cluster.organism_id)

    # Now go through the genotype dict and remove entries that only show up
    # once
    filtered_genotypes = {}
    for string_hash, item in genotypes.items():
        if item.count > min_occurence:
            filtered_genotypes[string_hash] = item

    # clear the <genotypes> dict
    genotypes = {}

    # replace it with a numbered version
    genotypes = dict(zip(range(0, len(filtered_genotypes)),
                         filtered_genotypes.values()))

    return genotypes


def generate_cytotable(genotypes):
    d = {}
    for key, item in genotypes.items():
        for accession in item.accessions:
            d[accession] = key

    return pd.DataFrame.from_dict(d, orient='index')


class analysis:
    def __init__(self, filefield, query):
        """
        Input:
        - a FileField object with a .filename property (string) that says
        what the file is called in the /uploads directory
        - a <query> object
        """
        self.query = query
        self.dataset = parse_dataset(importdata(filefield)[0],
                                     importdata(filefield)[1],
                                     query)
        self.genotypes = bin_hashes(self.dataset)
        self.cytotable = generate_cytotable(self.genotypes)

# TODO: format html report so its purdy
# TODO: implement hover/onclick info (accessions, gene names, organisms)
