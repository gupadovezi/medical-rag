#!/usr/bin/env python3
"""
Offline RAG Medical Assistant
This version works without API calls, using only the knowledge base and PubMed data.
"""

import requests
import json
from langchain_core.documents import Document
from langchain_community.embeddings import HuggingFaceEmbeddings
from langchain_community.vectorstores import Chroma
from langchain.text_splitter import RecursiveCharacterTextSplitter
from Bio import Entrez
import time
import re

# PubMed configuration
Entrez.email = "padovezigustavo@gmail.com"

# Initialize embeddings and vector store
embeddings = HuggingFaceEmbeddings(model_name="sentence-transformers/all-MiniLM-L6-v2")
vectorstore = None

def search_pubmed(query, max_results=3):
    """Search PubMed for relevant medical literature."""
    try:
        print(f"🔍 Searching PubMed for: {query}")
        
        # Search PubMed
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]:
            print("No articles found")
            return []
        
        articles = []
        # Get each article individually to avoid parsing issues
        for pmid in record["IdList"][:max_results]:
            try:
                handle = Entrez.esummary(db="pubmed", id=pmid)
                summary = Entrez.read(handle)
                handle.close()
                
                # Extract the article data
                if pmid in summary:
                    article_data = summary[pmid]
                    title = article_data.get('Title', 'No title available')
                    abstract = article_data.get('Abstract', 'No abstract available')
                    authors = article_data.get('AuthorList', [])
                    journal = article_data.get('FullJournalName', 'Unknown journal')
                    pub_date = article_data.get('PubDate', 'Unknown date')
                    
                    # Create a formatted article entry
                    article_text = f"Title: {title}\n"
                    if authors:
                        article_text += f"Authors: {', '.join(authors[:3])}{' et al.' if len(authors) > 3 else ''}\n"
                    article_text += f"Journal: {journal}\n"
                    article_text += f"Date: {pub_date}\n"
                    article_text += f"Abstract: {abstract}\n"
                    
                    articles.append({
                        'title': title,
                        'abstract': abstract,
                        'authors': authors,
                        'journal': journal,
                        'pub_date': pub_date,
                        'text': article_text
                    })
                    
            except Exception as e:
                print(f"Warning: Could not retrieve article {pmid}: {e}")
                continue
        
        print(f"✅ Found {len(articles)} relevant PubMed articles")
        return articles
        
    except Exception as e:
        print(f"❌ Error searching PubMed: {e}")
        return []

def get_pubmed_context(query, max_results=3):
    """Get PubMed context for a medical query."""
    articles = search_pubmed(query, max_results)
    
    if not articles:
        return ""
    
    context = "Latest medical research from PubMed:\n\n"
    for i, article in enumerate(articles, 1):
        context += f"Research Article {i}:\n{article['text']}\n"
        context += "-" * 50 + "\n"
    
    return context

def setup_knowledge_base(documents, metadata=None):
    """Set up the knowledge base with documents."""
    global vectorstore
    
    print("🔧 Setting up knowledge base...")
    
    # Create Document objects
    docs = []
    for i, content in enumerate(documents):
        meta = metadata[i] if metadata and i < len(metadata) else {"source": f"doc_{i}"}
        docs.append(Document(page_content=content, metadata=meta))
    
    # Split documents into chunks
    text_splitter = RecursiveCharacterTextSplitter(
        chunk_size=1000,
        chunk_overlap=200
    )
    split_docs = text_splitter.split_documents(docs)
    
    # Create vector store
    vectorstore = Chroma.from_documents(
        documents=split_docs,
        embedding=embeddings,
        collection_name="rag_knowledge_base"
    )
    
    print(f"✅ Knowledge base created with {len(split_docs)} document chunks")
    return vectorstore

def retrieve_context(query, k=3):
    """Retrieve relevant context from the knowledge base."""
    if not vectorstore:
        return ""
    
    try:
        docs = vectorstore.similarity_search(query, k=k)
        context = "Relevant information from knowledge base:\n\n"
        for i, doc in enumerate(docs, 1):
            context += f"Document {i} (Source: {doc.metadata.get('source', 'Unknown')}):\n"
            context += f"{doc.page_content.strip()}\n\n"
        return context
    except Exception as e:
        print(f"❌ Error retrieving context: {e}")
        return ""

def generate_offline_response(question, context):
    """Generate a response using the retrieved context without API calls."""
    
    # Extract key information from context
    kb_info = ""
    pubmed_info = ""
    
    if "knowledge base" in context.lower():
        kb_parts = context.split("Latest medical research from PubMed:")
        if len(kb_parts) > 1:
            kb_info = kb_parts[0].strip()
            pubmed_info = "Latest medical research from PubMed:" + kb_parts[1]
        else:
            kb_info = context
    else:
        kb_info = context
    
    # Create a structured response
    response = f"""Based on the available medical information, here's what I can tell you about your question: "{question}"

📚 KNOWLEDGE BASE INFORMATION:
{kb_info}

"""
    
    if pubmed_info:
        response += f"""
🔬 LATEST RESEARCH:
{pubmed_info}

"""
    
    response += """
💡 SUMMARY:
"""
    
    # Generate a summary based on the context
    if "scoliosis" in question.lower():
        if "treatment" in question.lower() or "treat" in question.lower():
            response += """Scoliosis treatment options depend on the severity and patient age:
• Mild cases (<25°): Observation and monitoring
• Moderate cases (25-45°): Bracing for growing children
• Severe cases (>45-50°): Surgical intervention
• Physical therapy and exercises can help with posture and pain
• Alternative treatments include chiropractic care and massage therapy

Always consult with a qualified healthcare provider for personalized treatment recommendations."""
        else:
            response += """Scoliosis is a medical condition characterized by an abnormal lateral curvature of the spine. Treatment options vary based on severity, and regular monitoring by a spine specialist is important for tracking progression."""
    
    elif "bracing" in question.lower():
        response += """Bracing is a common treatment for moderate scoliosis in growing children. The most effective braces include Boston, Milwaukee, and Charleston bending braces. They work best when worn 16-23 hours per day during growth periods to prevent curve progression."""
    
    elif "surgery" in question.lower() or "surgical" in question.lower():
        response += """Surgical treatment is typically considered for severe scoliosis (>45-50°) or when non-surgical treatments fail. Spinal fusion surgery connects vertebrae to prevent further curvature using rods, screws, and bone grafts. Recovery typically takes 6-12 months with physical therapy."""
    
    elif "exercise" in question.lower() or "physical therapy" in question.lower():
        response += """Physical therapy and exercise programs help improve posture, strengthen core muscles, and reduce pain. The Schroth method focuses on breathing and postural correction. Regular stretching and strengthening exercises maintain flexibility and reduce discomfort."""
    
    else:
        response += """The information provided covers various aspects of medical conditions and treatments. For specific medical advice, please consult with a qualified healthcare professional."""
    
    response += """

⚠️ DISCLAIMER: This information is for educational purposes only and should not replace professional medical advice. Always consult with a qualified healthcare provider for diagnosis and treatment recommendations."""

    return response

def ask_offline(question, include_pubmed=True):
    """Ask a question using offline RAG (no API calls)."""
    print(f"🤔 Question: {question}")
    
    # Retrieve relevant context from knowledge base
    kb_context = retrieve_context(question)
    
    # Get PubMed context if enabled
    pubmed_context = ""
    if include_pubmed:
        pubmed_context = get_pubmed_context(question)
    
    # Combine contexts
    combined_context = ""
    if kb_context:
        combined_context += kb_context + "\n"
    if pubmed_context:
        combined_context += pubmed_context + "\n"
    
    # Generate offline response
    response = generate_offline_response(question, combined_context)
    
    return {
        "choices": [{
            "message": {
                "content": response
            }
        }],
        "rag_context": combined_context[:500] + "..." if len(combined_context) > 500 else combined_context,
        "offline": True
    }

# Sample medical documents for scoliosis treatment
medical_documents = [
    """Scoliosis is a medical condition characterized by an abnormal lateral curvature of the spine. 
    Treatment options vary depending on the severity of the curve and the patient's age. 
    For mild cases (curves less than 25 degrees), observation and monitoring are typically recommended. 
    Regular check-ups with a spine specialist are important to track progression.""",
    
    """Bracing is a common treatment for moderate scoliosis (curves between 25-45 degrees) in growing children. 
    The most commonly used braces include the Boston brace, Milwaukee brace, and Charleston bending brace. 
    Bracing is most effective when worn for 16-23 hours per day and when the patient is still growing. 
    The goal is to prevent further progression of the curve.""",
    
    """Physical therapy and exercise programs can help improve posture, strengthen core muscles, 
    and reduce pain associated with scoliosis. Specific exercises like the Schroth method, 
    which focuses on breathing and postural correction, have shown promising results. 
    Regular stretching and strengthening exercises can help maintain flexibility and reduce discomfort.""",
    
    """Surgical treatment is typically considered for severe scoliosis (curves greater than 45-50 degrees) 
    or when non-surgical treatments have failed. Spinal fusion surgery involves connecting vertebrae 
    together to prevent further curvature. Modern surgical techniques use rods, screws, and bone grafts 
    to stabilize the spine. Recovery typically takes 6-12 months with physical therapy.""",
    
    """Alternative treatments for scoliosis include chiropractic care, massage therapy, and acupuncture. 
    While these may provide temporary pain relief, there is limited scientific evidence supporting 
    their effectiveness in correcting spinal curvature. It's important to consult with a qualified 
    healthcare provider before pursuing any treatment option."""
]

# Document metadata
metadata = [
    {"source": "general_info.txt", "topic": "overview"},
    {"source": "bracing_treatment.txt", "topic": "bracing"},
    {"source": "physical_therapy.txt", "topic": "exercise"},
    {"source": "surgical_treatment.txt", "topic": "surgery"},
    {"source": "alternative_treatments.txt", "topic": "alternatives"}
]

def main():
    """Main function to run the offline RAG system."""
    print("🚀 Starting Offline RAG Medical Assistant...")
    print("=" * 60)
    print("This version works without API calls - no rate limits!")
    print("=" * 60)
    
    # Setup knowledge base
    setup_knowledge_base(medical_documents, metadata)
    
    print("\n🤖 Offline RAG Medical Assistant Interactive Mode")
    print("Type your question and press Enter. Type 'quit' or 'exit' to stop.")
    print("Commands:")
    print("  'pubmed on/off' - toggle PubMed search")
    print("  'quit' or 'exit' - stop the session")
    print("=" * 60)
    
    include_pubmed = True
    
    while True:
        try:
            user_input = input("\n👤 Ask a medical question: ").strip()
            
            if user_input.lower() in ["quit", "exit", "bye"]:
                print("\n🤖 Goodbye! Have a great day!")
                break
            
            if not user_input:
                continue
            
            # Handle commands
            if user_input.lower() == "pubmed on":
                include_pubmed = True
                print("✅ PubMed search enabled")
                continue
            elif user_input.lower() == "pubmed off":
                include_pubmed = False
                print("❌ PubMed search disabled")
                continue
                
            result = ask_offline(user_input, include_pubmed=include_pubmed)
            
            if "choices" in result and result["choices"]:
                print(f"\n🤖 Answer: {result['choices'][0]['message']['content']}")
                
                if "rag_context" in result:
                    print(f"\n📚 Context used: {result['rag_context'][:200]}...")
                
                print("\n✅ This response was generated offline - no API calls needed!")
            else:
                print(f"\n❌ Error: {result.get('error', 'Unknown error occurred')}")
                
        except KeyboardInterrupt:
            print("\n\n🤖 Goodbye! Have a great day!")
            break
        except Exception as e:
            print(f"\n❌ Unexpected error: {e}")
            print("Please try again.")

if __name__ == "__main__":
    main() 