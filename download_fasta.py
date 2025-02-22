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
        search_handle = Entrez.esearch(db="nucleotide", term=f"{species_name}[ORGN] AND complete genome", retmax=1)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        if not search_results["IdList"]:
            print("No reference genome found for the given species.")
            return None

        accession = search_results["IdList"][0]
        print(f"Found Reference Genome Accession: {accession}")
        return accession

    except Exception as e:
        print(f"Error during species search: {e}")
        return None

def download_fasta(accession, output_file):
    """
    Downloads a FASTA file from NCBI using an accession number.
    
    Parameters:
        accession (str): The NCBI GenBank accession number.
        output_file (str): The name of the output file to save the FASTA data.
    """
    try:
        print(f"Fetching FASTA record for: {accession}")
        
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()

        with open(output_file, "w") as file:
            file.write(fasta_data)
        
        print(f"Download complete: {output_file}")

    except Exception as e:
        print(f"Error fetching FASTA file: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download reference FASTA files from NCBI using a species name or GenBank accession ID.")
    parser.add_argument("-i", "--input", required=True, help="Species name or GenBank accession ID")

    args = parser.parse_args()
    input_query = args.input

    if is_accession(input_query):
        # Input is an accession ID, download FASTA directly
        output_filename = f"{input_query}.fasta"
        download_fasta(input_query, output_filename)
    else:
        # Input is a species name, search for the reference genome
        accession_number = search_species(input_query)
        if accession_number:
            output_filename = f"{input_query.replace(' ', '_')}.fasta"
            download_fasta(accession_number, output_filename)