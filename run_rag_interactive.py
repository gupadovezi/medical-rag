#!/usr/bin/env python3
"""
Interactive RAG Medical Assistant with PubMed Integration
Run this script in a new terminal window to start the interactive RAG system.
"""

import requests
import json
from langchain_core.documents import Document
from langchain_community.embeddings import HuggingFaceEmbeddings
from langchain_community.vectorstores import Chroma
from langchain.text_splitter import RecursiveCharacterTextSplitter
from Bio import Entrez
import time

# Your API configuration
API_KEY = "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"
MODEL_ID = "deepseek/deepseek-v3-base:free"  # Default model
ALTERNATIVE_MODELS = [
    "anthropic/claude-3-haiku:free",
    "meta-llama/llama-3.1-8b-instruct:free",
    "google/gemini-flash-1.5:free"
]
current_model_index = 0

# PubMed configuration (required by NCBI)
Entrez.email = "padovezigustavo@gmail.com"  # Replace with your email for PubMed access

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

def ask_with_rag(question, temperature=0.7, max_tokens=1024, include_pubmed=True):
    """Ask a question with RAG enhancement and PubMed integration."""
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
    
    # Create RAG-enhanced prompt
    if combined_context:
        enhanced_prompt = f"""You are a helpful AI assistant with access to a knowledge base and the latest medical research from PubMed. 
Use the following relevant information to answer the user's question accurately and comprehensively.

{combined_context}

User Question: {question}

Please provide a detailed answer based on the information above. If the information doesn't fully address the question, acknowledge what you know and what might need additional research. When citing research, mention the source (knowledge base or PubMed research)."""
    else:
        enhanced_prompt = question
    
    # Prepare messages
    messages = [
        {"role": "system", "content": "You are a helpful AI assistant with access to medical knowledge and research."},
        {"role": "user", "content": enhanced_prompt}
    ]
    
    # Make API call with retry logic
    max_retries = 3
    base_delay = 5  # Increased delay
    
    for attempt in range(max_retries):
        try:
            headers = {
                "Authorization": f"Bearer {API_KEY}",
                "HTTP-Referer": f"https://openrouter.ai/{MODEL_ID}",
                "X-Title": "RAG-Enhanced LLM",
                "Content-Type": "application/json"
            }
            
            payload = {
                "model": MODEL_ID,
                "messages": messages,
                "temperature": temperature,
                "max_tokens": max_tokens,
                "stream": False,
                "top_p": 0.95
            }
            
            response = requests.post(
                "https://openrouter.ai/api/v1/chat/completions",
                headers=headers,
                data=json.dumps(payload)
            )
            
            if response.status_code == 429:
                delay = base_delay * (2 ** attempt)  # Exponential backoff
                print(f"⚠️ Rate limit hit. Waiting {delay} seconds before retry {attempt + 1}/{max_retries}...")
                time.sleep(delay)
                continue
            elif response.status_code != 200:
                response.raise_for_status()
            
            result = response.json()
            
            # Add context information
            if combined_context and "choices" in result:
                result["rag_context"] = combined_context[:500] + "..." if len(combined_context) > 500 else combined_context
            
            return result
            
        except requests.exceptions.RequestException as e:
            if "429" in str(e) and attempt < max_retries - 1:
                delay = base_delay * (2 ** attempt)
                print(f"⚠️ Rate limit error. Waiting {delay} seconds before retry {attempt + 1}/{max_retries}...")
                time.sleep(delay)
                continue
            else:
                print(f"❌ API call failed: {e}")
                return {"error": str(e)}
    
    # If all retries failed, provide a fallback response
    print("⚠️ All retry attempts failed due to rate limiting. Providing fallback response...")
    return provide_fallback_response(question, combined_context)

def provide_fallback_response(question, context):
    """Provide a fallback response when API is rate limited."""
    fallback_response = f"""I apologize, but I'm currently experiencing high demand and cannot access the AI model right now. 

However, I can provide you with the relevant information I found:

{context}

Based on this information, here's what I can tell you about your question: "{question}"

For a more detailed AI-generated response, please try again in a few minutes when the service is less busy.

You can also:
- Try asking a different question
- Use 'pubmed off' to disable PubMed search for faster responses
- Wait a few minutes and try again"""
    
    return {
        "choices": [{
            "message": {
                "content": fallback_response
            }
        }],
        "rag_context": context[:500] + "..." if len(context) > 500 else context,
        "fallback": True
    }

def switch_model():
    """Switch to the next available model."""
    global current_model_index, MODEL_ID
    current_model_index = (current_model_index + 1) % len(ALTERNATIVE_MODELS)
    MODEL_ID = ALTERNATIVE_MODELS[current_model_index]
    print(f"✅ Switched to model: {MODEL_ID}")

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
    """Main function to run the interactive RAG system."""
    print("🚀 Starting RAG Medical Assistant with PubMed Integration...")
    print("=" * 60)
    
    # Setup knowledge base
    setup_knowledge_base(medical_documents, metadata)
    
    print("\n🤖 RAG Medical Assistant Interactive Mode")
    print("Type your question and press Enter. Type 'quit' or 'exit' to stop.")
    print("Commands:")
    print("  'pubmed on/off' - toggle PubMed search")
    print("  'switch model' - try a different AI model")
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
            elif user_input.lower() == "switch model":
                switch_model()
                continue
                
            result = ask_with_rag(user_input, include_pubmed=include_pubmed)
            
            if "choices" in result and result["choices"]:
                print(f"\n🤖 Answer: {result['choices'][0]['message']['content']}")
                
                if "rag_context" in result:
                    print(f"\n📚 Context used: {result['rag_context'][:200]}...")
                
                if result.get("fallback"):
                    print("\n⚠️ This was a fallback response due to rate limiting.")
            else:
                print(f"\n❌ Error: {result.get('error', 'Unknown error occurred')}")
            
            # Add a small delay to prevent rate limiting
            time.sleep(1)
                
        except KeyboardInterrupt:
            print("\n\n🤖 Goodbye! Have a great day!")
            break
        except Exception as e:
            print(f"\n❌ Unexpected error: {e}")
            print("Please try again.")

if __name__ == "__main__":
    main() 