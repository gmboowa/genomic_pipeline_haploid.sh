import ssl
import argparse
from Bio import Entrez

# Fix SSL issue for NCBI access
ssl._create_default_https_context = ssl._create_unverified_context

# Set your email for NCBI access (replace with your email)
Entrez.email = "your_email@example.com"

def is_accession(input_str):
    """
    Checks if the input is a GenBank accession ID.
    A valid accession usually contains letters followed by numbers (e.g., KR781608.1).
    """
    return any(char.isdigit() for char in input_str) and input_str[0].isalpha()

def search_species(species_name):
    """
    Searches for a species in the NCBI Nucleotide database and retrieves the first matching GenBank accession.
    
    Parameters:
        species_name (str): The scientific name of the species.
    
    Returns:
        str: The first GenBank accession number found.
    """
    try:
        print(f"Searching for species: {species_name}")
        search_handle = Entrez.esearch(db="nucleotide", term=f"{species_name}[ORGN]", retmax=1)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        if not search_results["IdList"]:
            print("No GenBank records found for the given species.")
            return None

        accession = search_results["IdList"][0]
        print(f"Found GenBank Accession: {accession}")
        return accession

    except Exception as e:
        print(f"Error during species search: {e}")
        return None

def download_genbank(accession, output_file):
    """
    Downloads a GenBank file from NCBI using an accession number.
    
    Parameters:
        accession (str): The NCBI GenBank accession number.
        output_file (str): The name of the output file to save the GenBank data.
    """
    try:
        print(f"Fetching GenBank record for: {accession}")
        
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        genbank_data = handle.read()
        handle.close()

        with open(output_file, "w") as file:
            file.write(genbank_data)
        
        print(f"Download complete: {output_file}")

    except Exception as e:
        print(f"Error fetching GenBank file: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download GenBank files from NCBI using a species name or accession ID.")
    parser.add_argument("-i", "--input", required=True, help="Species name or GenBank accession ID")

    args = parser.parse_args()
    input_query = args.input

    if is_accession(input_query):
        # Input is an accession ID, download directly
        output_filename = f"{input_query}.gb"
        download_genbank(input_query, output_filename)
    else:
        # Input is a species name, search for the first available record
        accession_number = search_species(input_query)
        if accession_number:
            output_filename = f"{input_query.replace(' ', '_')}.gb"
            download_genbank(accession_number, output_filename)
