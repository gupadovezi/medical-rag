#!/usr/bin/env python3
"""
Test script for PubMed integration
"""

from Bio import Entrez

# Set email
Entrez.email = "padovezigustavo@gmail.com"

def test_pubmed_search():
    """Test PubMed search functionality."""
    try:
        print("🔍 Testing PubMed search...")
        
        # Simple search
        handle = Entrez.esearch(db="pubmed", term="scoliosis treatment", retmax=2)
        record = Entrez.read(handle)
        handle.close()
        
        print(f"✅ Search successful! Found {len(record['IdList'])} articles")
        
        if record["IdList"]:
            # Test individual article retrieval
            pmid = record["IdList"][0]
            handle = Entrez.esummary(db="pubmed", id=pmid)
            summary = Entrez.read(handle)
            handle.close()
            
            if pmid in summary:
                article_data = summary[pmid]
                print(f"✅ Article retrieval successful!")
                print(f"Title: {article_data.get('Title', 'No title')}")
                print(f"Journal: {article_data.get('FullJournalName', 'Unknown')}")
            
        return True
        
    except Exception as e:
        print(f"❌ Error: {e}")
        return False

if __name__ == "__main__":
    test_pubmed_search() 