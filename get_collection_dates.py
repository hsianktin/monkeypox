# Import BioPython modules
from Bio import Entrez
from Bio import SeqIO
from datetime import datetime
import csv
import pandas as pd
import re

def convert_date(date_str):
  # Try to match different date formats using regular expressions
  match1 = re.match(r'(\d{4})-(\d{2})-(\d{2})', date_str) # YYYY-MM-DD
  match2 = re.match(r'(\d{4})-(\d{2})', date_str) # YYYY-MM
  match3 = re.match(r'(\d{4})', date_str) # YYYY
  match4 = re.match(r'(\d{2})-(\w{3})-(\d{4})', date_str) # DD-MMM-YYYY
  match5 = re.match(r'(\w{3})-(\d{4})', date_str) # MMM-YYYY

  if match1:
    return date_str # Already in desired format
  elif match2:
    return f'{date_str}-01' # Add day as 01
  elif match3:
    return f'{date_str}-01-01' # Add month and day as 01
  elif match4:
    # Convert month name to number using datetime.strptime and strftime methods[^1^][1] [^2^][2]
    month_num = datetime.strptime(match4.group(2), '%b').strftime('%m')
    return f'{match4.group(3)}-{month_num}-{match4.group(1)}' # Rearrange year, month and day
  elif match5:
    # Convert month name to number using datetime.strptime and strftime methods[^1^][1] [^2^][2]
    month_num = datetime.strptime(match5.group(1), '%b').strftime('%m')
    return f'{match5.group(2)}-{month_num}-01' # Add day as 01

# Set your email address for NCBI
Entrez.email = "xiangting.li@ucla.edu"

# read metadata.csv file and get a list of accession numbers
df = pd.read_csv('accessions.csv')
accessions = df['Accession'].to_list()

# read meta_extra_data.csv file and get a list of accession numbers
df = pd.read_csv('accessions_extra.csv')
extra_accessions = df['Accession'].to_list()

# remove duplicates of the form 'MT123456.1' and 'MT123456'
# in the extra_accessions list, the accession numbers are in the form 'MT123456.1'
# in the accessions list, the accession numbers are in the form 'MT123456'
# so we remove the duplicates
extra_accessions = [acc[:-2] for acc in extra_accessions if acc[:-2] not in accessions]

# combine the two lists
accessions = accessions + extra_accessions



# Create an empty list to store accs with collection dates 
records_with_dates = []

# Define a function to get collection date from accession number
def get_collection_date(accession):
    # Use Entrez.efetch to fetch a GenBank record by accession number
    handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
    # Parse the record using SeqIO
    record = SeqIO.read(handle, "genbank")
    # Close the handle
    handle.close()
    # Loop through the features of the record
    for feature in record.features:
        # If the feature is a source feature
        if feature.type == "source":
            # Try to get the collection_date qualifier
            try:
                collection_date = feature.qualifiers["collection_date"][0]
                collection_date = convert_date(collection_date)
                # Return the date in yyyy-mm-dd format if possible
                try:
                    collection_date_dt = datetime.strptime(collection_date, "%Y-%m-%d")
                except:
                    try:
                        collection_date_dt =  datetime.strptime(collection_date, "%b-%Y")
                    except:
                        collection_date_dt = datetime.strptime(collection_date, "%Y")

                # Convert datetime object to yyyy-mm-dd format using strftime method 
                collection_date_ymd = collection_date_dt.strftime("%Y-%m-%d")
                # Print collection date in yyyy-mm-dd format 
                # print(f"Collection date: {collection_date_ymd}")
                return collection_date_ymd
            except KeyError:
                # If there is no collection_date qualifier, return NAN
                return "NAN"

# for acc in accessions:
#     date = get_collection_date(acc)
#     if date != "NAN":
#         # store both acc and date in a list
#         records_with_dates.append([acc, date])

# this take some time, so we use progressbar of tqdm
from tqdm import tqdm
for acc in tqdm(accessions):
    date = get_collection_date(acc)
    if date != "NAN":
        # store both acc and date in a list
        records_with_dates.append([acc, date])

# Write the list to a csv file by converting it to a dataframe
df = pd.DataFrame(records_with_dates, columns=['Accession', 'Date'])
df.to_csv('collection_dates.csv', index=False)