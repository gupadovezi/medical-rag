import math
import os
from pathlib import Path
import PyPDF2
import pandas as pd
from datetime import datetime
import re
from typing import List, Dict, Any

# Basic math operations
print(f"Pi: {math.pi}")
print(f"Square root of 16: {math.sqrt(16)}")
print(f"Cosine of 0: {math.cos(0)}")
print(f"Factorial of 5: {math.factorial(5)}")

def read_protocol(protocol_path):
    """Read the systematic review protocol to understand extraction requirements."""
    try:
        with open(protocol_path, 'rb') as file:
            reader = PyPDF2.PdfReader(file)
            protocol_text = ""
            for page in reader.pages:
                protocol_text += page.extract_text() + "\n"
            return protocol_text
    except Exception as e:
        print(f"Error reading protocol: {str(e)}")
        return ""

def extract_study_info(text, filename):
    """Extract specific study information from the text based on systematic review protocol."""
    # Initialize dictionary with default values
    study_info = {
        'Study ID': '',
        'Reference': '',
        'Study Design': '',
        'Year': '',
        'Language': '',
        'Country': '',
        'Clinical Trial Registration': '',
        'Sample Size': '',
        'Population': '',
        'Intervention': '',
        'Comparator': '',
        'Outcomes': '',
        'Effect Size': '',
        'Measures of Treatment Effect': '',
        'Duration of Participation': '',
        'Missing Data': '',
        'Reason for Missing Data': '',
        'Inclusion Criteria': '',
        'Baseline Pain': '',
        'Age (Mean)': '',
        'Age (SD)': '',
        'Sex (M/F)': '',
        'Pain Duration': '',
        'Intervention Frequency': '',
        'Intervention Duration': '',
        'Follow-up Duration': '',
        'Outcome Assessment Timepoints': '',
        'Risk of Bias': '',
        'Funding Source': '',
        'Conflicts of Interest': ''
    }
    
    # Extract year from filename or text
    year_match = re.search(r'20\d{2}', filename)
    if year_match:
        study_info['Year'] = year_match.group(0)
    
    # Extract study ID from filename
    study_id_match = re.search(r'(\d+)\s*-', filename)
    if study_id_match:
        study_info['Study ID'] = study_id_match.group(1)
    
    # Extract study design
    design_patterns = [
        r'(randomized controlled trial|RCT|randomised controlled trial)',
        r'(cluster randomized|cluster randomised)',
        r'(crossover|cross-over)',
        r'(non-randomized|non-randomised)',
        r'(quasi-experimental)',
        r'(pilot study|pilot trial)'
    ]
    for pattern in design_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            study_info['Study Design'] = match.group(1)
            break
    
    # Extract sample size
    sample_size_patterns = [
        r'n\s*=\s*(\d+)',
        r'sample size.*?(\d+)',
        r'(\d+)\s*participants',
        r'(\d+)\s*patients',
        r'(\d+)\s*subjects'
    ]
    for pattern in sample_size_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            study_info['Sample Size'] = match.group(1)
            break
    
    # Extract age information
    age_patterns = [
        r'mean age.*?(\d+(?:\.\d+)?)\s*(?:years|y\.o\.|y\.o)',
        r'aged.*?(\d+(?:\.\d+)?)\s*(?:years|y\.o\.|y\.o)',
        r'age.*?(\d+(?:\.\d+)?)\s*(?:years|y\.o\.|y\.o)'
    ]
    for pattern in age_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            study_info['Age (Mean)'] = match.group(1)
            break
    
    # Extract age SD
    age_sd_patterns = [
        r'age.*?SD\s*=\s*(\d+(?:\.\d+)?)',
        r'age.*?standard deviation.*?(\d+(?:\.\d+)?)',
        r'age.*?\(\s*(\d+(?:\.\d+)?)\s*\)'
    ]
    for pattern in age_sd_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            study_info['Age (SD)'] = match.group(1)
            break
    
    # Extract sex ratio
    sex_patterns = [
        r'(\d+)\s*male.*?(\d+)\s*female',
        r'(\d+)\s*M.*?(\d+)\s*F',
        r'(\d+)\s*men.*?(\d+)\s*women'
    ]
    for pattern in sex_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            study_info['Sex (M/F)'] = f"{match.group(1)}/{match.group(2)}"
            break
    
    # Extract intervention details
    intervention_patterns = [
        r'intervention.*?([^.]*?\.)',
        r'treatment.*?([^.]*?\.)',
        r'therapy.*?([^.]*?\.)'
    ]
    for pattern in intervention_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            study_info['Intervention'] = match.group(1).strip()
            break
    
    # Extract comparator
    comparator_patterns = [
        r'control group.*?([^.]*?\.)',
        r'comparison group.*?([^.]*?\.)',
        r'versus.*?([^.]*?\.)',
        r'compared to.*?([^.]*?\.)'
    ]
    for pattern in comparator_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            study_info['Comparator'] = match.group(1).strip()
            break
    
    # Extract outcomes
    outcome_patterns = [
        r'primary outcome.*?([^.]*?\.)',
        r'secondary outcome.*?([^.]*?\.)',
        r'outcomes.*?([^.]*?\.)'
    ]
    for pattern in outcome_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            study_info['Outcomes'] = match.group(1).strip()
            break
    
    # Extract inclusion criteria
    inclusion_patterns = [
        r'inclusion criteria.*?([^.]*?\.)',
        r'eligible.*?([^.]*?\.)',
        r'included.*?([^.]*?\.)'
    ]
    for pattern in inclusion_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            study_info['Inclusion Criteria'] = match.group(1).strip()
            break
    
    # Extract baseline pain
    pain_patterns = [
        r'baseline pain.*?(\d+(?:\.\d+)?)',
        r'initial pain.*?(\d+(?:\.\d+)?)',
        r'pain score.*?(\d+(?:\.\d+)?)'
    ]
    for pattern in pain_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            study_info['Baseline Pain'] = match.group(1)
            break
    
    # Extract duration of pain
    pain_duration_patterns = [
        r'pain duration.*?(\d+(?:\.\d+)?)\s*(?:years|months|weeks)',
        r'chronic pain.*?(\d+(?:\.\d+)?)\s*(?:years|months|weeks)',
        r'pain for.*?(\d+(?:\.\d+)?)\s*(?:years|months|weeks)'
    ]
    for pattern in pain_duration_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            study_info['Pain Duration'] = match.group(1)
            break
    
    # Extract intervention duration
    intervention_duration_patterns = [
        r'intervention.*?(\d+(?:\.\d+)?)\s*(?:weeks|months|sessions)',
        r'treatment.*?(\d+(?:\.\d+)?)\s*(?:weeks|months|sessions)',
        r'therapy.*?(\d+(?:\.\d+)?)\s*(?:weeks|months|sessions)'
    ]
    for pattern in intervention_duration_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            study_info['Intervention Duration'] = match.group(1)
            break
    
    # Extract follow-up duration
    followup_patterns = [
        r'follow-up.*?(\d+(?:\.\d+)?)\s*(?:weeks|months|years)',
        r'follow up.*?(\d+(?:\.\d+)?)\s*(?:weeks|months|years)',
        r'followed.*?(\d+(?:\.\d+)?)\s*(?:weeks|months|years)'
    ]
    for pattern in followup_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            study_info['Follow-up Duration'] = match.group(1)
            break
    
    return study_info

def read_pdf_file(pdf_path):
    """Read and extract text from a PDF file."""
    try:
        with open(pdf_path, 'rb') as file:
            reader = PyPDF2.PdfReader(file)
            text = ""
            for page in reader.pages:
                text += page.extract_text() + "\n"
            
            # Extract study information
            study_info = extract_study_info(text, pdf_path.name)
            study_info['filename'] = pdf_path.name
            study_info['full_text'] = text
            
            return study_info
    except Exception as e:
        return {
            'filename': pdf_path.name,
            'error': f"Error reading file: {str(e)}"
        }

def process_pdf_directory(directory_path: str, protocol_path: str = None) -> List[Dict[str, Any]]:
    """
    Process all PDF files in a directory and return extracted data.
    
    Args:
        directory_path (str): Path to the directory containing PDF files
        protocol_path (str, optional): Path to the protocol file. If None, will use default extraction.
    
    Returns:
        List[Dict[str, Any]]: List of dictionaries containing extracted data
    """
    results = []
    pdf_files = list(Path(directory_path).glob("*.pdf"))
    
    for pdf_file in pdf_files:
        try:
            # Skip protocol file if it's in the same directory
            if protocol_path and str(pdf_file) == protocol_path:
                continue
                
            # Extract text from PDF
            text = extract_text_from_pdf(str(pdf_file))
            
            # Create result dictionary
            result = {
                'filename': pdf_file.name,
                'text': text,
                'processed_at': datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
            
            results.append(result)
            
        except Exception as e:
            print(f"Error processing {pdf_file.name}: {str(e)}")
            continue
    
    return results

if __name__ == "__main__":
    pdf_directory = "/Users/gustavopadovezi/Desktop/pdfs"
    protocol_path = "/Users/gustavopadovezi/Library/Containers/net.whatsapp.WhatsApp/Data/tmp/documents/0D215BB4-CD5D-4F34-B142-C18DA0560EE7/SR_InterventionsChronicMSKPain_OSF.pdf"
    process_pdf_directory(pdf_directory, protocol_path) 