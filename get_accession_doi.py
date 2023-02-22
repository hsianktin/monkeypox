#!/usr/bin/env python3
from Bio import Entrez
from Bio import SeqIO
from datetime import datetime
import csv
import pandas as pd
import re
from title2bib.crossref import get_from_title
from doi2bib.crossref import get_bib
import bibtexparser
import tqdm
import Levenshtein

# Set your email address for NCBI
Entrez.email = "xiangting.li@ucla.edu"

# # sample accession number
# accession = 'MK783028'

# # record instance
# handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
# # Parse the record using SeqIO
# record = SeqIO.read(handle, "genbank")
# # Close the handle
# handle.close()

# the title to doi function is extremely unreliable

# title = record.annotations['references'][0].title

def get_doi_from_title(title):
    found = False
    doi = ""
    found, item = get_from_title(title, True)
    if found:
        if Levenshtein.ratio(title, item["title"]) > 0.9:
            if item["is_crossref"]:
                if "DOI" in item:
                    doi = item["DOI"]
        else:
            found = False
    return found, doi

def fname(bib_lib_entry):
    info0 = bib_lib_entry['ID'].split('_')
    text = ''
    for i in range(0,len(info0)):
        text = text + info0[i]
    # if 'month' in bib_lib_entry:
    #     text = text + bib_lib_entry['month']
    return text

def get_doi_cite_entry(title):
    found, doi = get_doi_from_title(title)
    bib_cite_entry = ""
    if found:
        bibtex = get_bib(doi)
        if bibtex[0]:
            bib_lib = bibtexparser.loads(bibtex[1])
            bib_lib_entry = bib_lib.entries[0]
            bib_cite_entry = fname(bib_lib_entry)
    return found, doi, bib_cite_entry

# test get_doi_cite_entry
# get_doi_cite_entry(title)
# (True, '10.3390/v7042168', 'Nakazawa2015apr')
def split_date(date_str):
  # Try to match different date formats using regular expressions
  # match the types YYYY, MMM-YYYY, DD-MMM-YYYY,
  # for the first case, return the date as is
  # for the second case, return MMM, YYYY respectively
    # for the third case, return MMM, YYYY respectively
    match1 = re.match(r'(\d{4})', date_str) # YYYY
    match2 = re.match(r'(\w{3})-(\d{4})', date_str) # MMM-YYYY
    match3 = re.match(r'(\d{2})-(\w{3})-(\d{4})', date_str) # DD-MMM-YYYY
    if match1:
        return date_str, date_str
    elif match2:
        return match2.group(1), match2.group(2)
    elif match3:
        return match3.group(2), match3.group(3)

# split_date('2020')
# ('2020', '2020')
# split_date('Apr-2020')
# ('Apr', '2020')
# split_date('01-Apr-2020')
# ('Apr', '2020')

def get_ref_from_accession(accession):
    found = False
    doi = ""
    bib_cite_entry = ""
    handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    try:
        title = record.annotations['references'][0].title
        found, doi, bib_cite_entry = get_doi_cite_entry(title)
    except:
        pass
    return found, doi, bib_cite_entry

def build_bibtex(accession):
    handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    try:
        title = record.annotations['references'][0].title
        journal = record.annotations['references'][0].journal
        authors = record.annotations['references'][0].authors
        fauthor = authors.split(',')[0]
        # replace space in fauthor with dash
        fauthor = fauthor.replace(' ', '-')
        date = record.annotations['date']
        month, year = split_date(date)
        id = fauthor + year
        entry = {
            'ID': id,
            'title': title,
            'journal': journal,
            'author': authors,
            'year': year,
            'ENTRYTYPE': 'article'
        }
        return entry, id
    except:
        pass
    return {}, ""

# test build_bibtex
# build_bibtex(accession)


# df = pd.read_csv('accessions.csv')
df = pd.read_csv('meta_data.csv')
accessions = df['Accession'].to_list()

# empty list to store dois and bib_cite_entries
dois = []
bib_cite_entries = []
bibtex_data = bibtexparser.bibdatabase.BibDatabase()


# loop through the accessions
# use tqdm to show the progress bar
for acc in tqdm.tqdm(accessions):
    found, doi, bib_cite_entry = get_ref_from_accession(acc)
    bibtex_data_entry, id = build_bibtex(acc)
    bibtex_data.entries.append(bibtex_data_entry)
    dois.append(doi)
    if found:
        bib_cite_entries.append(bib_cite_entry)
    else:
        bib_cite_entries.append(id)

# add dois and bib_cite_entries to the dataframe
df['DOI'] = dois
df['Bib_cite_entry'] = bib_cite_entries

# save the dataframe to a csv file
df.to_csv('meta_data_doi.csv', index=False)

# remove duplicates in the bibtex data
bibtex_data.entries = list({v['ID']:v for v in bibtex_data.entries}.values())
# save the bibtex data to a bibtex file
with open('bibtex_data.bib', 'w') as bibtex_file:
    bibtexparser.dump(bibtex_data, bibtex_file)