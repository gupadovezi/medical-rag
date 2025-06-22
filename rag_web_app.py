import streamlit as st
import os
import json
import time
from typing import List, Dict, Any
import requests
from datetime import datetime
import chromadb
from chromadb.config import Settings

# Updated imports for the new langchain-huggingface package
from langchain_huggingface import HuggingFaceEmbeddings
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.vectorstores import Chroma
from langchain_core.documents import Document
from langchain.retrievers import ContextualCompressionRetriever
from langchain.retrievers.document_compressors import LLMChainExtractor
from Bio import Entrez
import re

# Page configuration
st.set_page_config(
    page_title="Medical AI Assistant",
    page_icon="🏥",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main-header {
        font-size: 3rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .sub-header {
        font-size: 1.5rem;
        color: #ff7f0e;
        margin-bottom: 1rem;
    }
    .info-box {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #1f77b4;
    }
    .response-box {
        background-color: #e8f4fd;
        padding: 1.5rem;
        border-radius: 0.5rem;
        border-left: 4px solid #ff7f0e;
        margin: 1rem 0;
    }
    .pubmed-box {
        background-color: #f0f8ff;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #2ca02c;
        margin: 0.5rem 0;
    }
</style>
""", unsafe_allow_html=True)

# Initialize session state
if 'vectorstore' not in st.session_state:
    st.session_state.vectorstore = None
if 'embeddings' not in st.session_state:
    st.session_state.embeddings = None
if 'chat_history' not in st.session_state:
    st.session_state.chat_history = []

# PubMed configuration
Entrez.email = "padovezigustavo@gmail.com"

def initialize_embeddings():
    """Initialize embeddings if not already done."""
    if st.session_state.embeddings is None:
        with st.spinner("Loading AI models..."):
            st.session_state.embeddings = HuggingFaceEmbeddings(
                model_name="sentence-transformers/all-MiniLM-L6-v2"
            )

def search_pubmed(query, max_results=3):
    """Search PubMed for relevant medical literature."""
    try:
        with st.spinner(f"🔍 Searching PubMed for: {query}"):
            # Search PubMed
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            
            if not record["IdList"]:
                return []
            
            articles = []
            # Get each article individually
            for pmid in record["IdList"][:max_results]:
                try:
                    handle = Entrez.esummary(db="pubmed", id=pmid)
                    summary = Entrez.read(handle)
                    handle.close()
                    
                    if pmid in summary:
                        article_data = summary[pmid]
                        title = article_data.get('Title', 'No title available')
                        abstract = article_data.get('Abstract', 'No abstract available')
                        authors = article_data.get('AuthorList', [])
                        journal = article_data.get('FullJournalName', 'Unknown journal')
                        pub_date = article_data.get('PubDate', 'Unknown date')
                        
                        articles.append({
                            'title': title,
                            'abstract': abstract,
                            'authors': authors,
                            'journal': journal,
                            'pub_date': pub_date
                        })
                        
                except Exception as e:
                    continue
            
            return articles
            
    except Exception as e:
        st.error(f"Error searching PubMed: {e}")
        return []

def setup_knowledge_base():
    """Set up the knowledge base with comprehensive medical documents."""
    if st.session_state.vectorstore is not None:
        return st.session_state.vectorstore
    
    # Comprehensive medical documents for musculoskeletal conditions
    medical_documents = [
        # Scoliosis
        Document(
            page_content="""Scoliosis is a medical condition characterized by an abnormal lateral curvature of the spine. 
            Treatment options vary depending on the severity of the curve and the patient's age. 
            For mild cases (curves less than 25 degrees), observation and monitoring are typically recommended. 
            Regular check-ups with a spine specialist are important to track progression.""",
            metadata={"source": "scoliosis_general.txt", "topic": "scoliosis", "category": "spine"}
        ),
        
        # Back Pain
        Document(
            page_content="""Lower back pain (lumbago) is one of the most common musculoskeletal conditions affecting adults. 
            Causes include muscle strain, herniated discs, spinal stenosis, and degenerative disc disease. 
            Treatment options include physical therapy, pain management, exercise, and in severe cases, surgery. 
            Prevention through proper posture, regular exercise, and ergonomic workplace setup is crucial.""",
            metadata={"source": "back_pain.txt", "topic": "back_pain", "category": "spine"}
        ),
        
        # Neck Pain
        Document(
            page_content="""Neck pain (cervicalgia) can result from poor posture, muscle strain, whiplash, or degenerative conditions. 
            Common causes include cervical spondylosis, herniated discs, and muscle tension. 
            Treatment includes physical therapy, neck exercises, pain management, and ergonomic adjustments. 
            Severe cases may require cervical spine surgery.""",
            metadata={"source": "neck_pain.txt", "topic": "neck_pain", "category": "spine"}
        ),
        
        # Arthritis
        Document(
            page_content="""Arthritis is inflammation of one or more joints, causing pain and stiffness. 
            Osteoarthritis is the most common type, affecting cartilage breakdown. 
            Rheumatoid arthritis is an autoimmune condition affecting joint linings. 
            Treatment includes medication, physical therapy, lifestyle changes, and joint replacement surgery for severe cases.""",
            metadata={"source": "arthritis.txt", "topic": "arthritis", "category": "joints"}
        ),
        
        # Sciatica
        Document(
            page_content="""Sciatica is pain that radiates along the path of the sciatic nerve, from the lower back down the leg. 
            It's usually caused by a herniated disc or bone spur compressing the nerve. 
            Symptoms include sharp pain, numbness, and tingling in the leg. 
            Treatment includes physical therapy, medication, epidural injections, and surgery for severe cases.""",
            metadata={"source": "sciatica.txt", "topic": "sciatica", "category": "nerves"}
        ),
        
        # Carpal Tunnel Syndrome
        Document(
            page_content="""Carpal tunnel syndrome is compression of the median nerve in the wrist, causing pain, numbness, and tingling in the hand. 
            Risk factors include repetitive hand movements, pregnancy, and certain medical conditions. 
            Treatment includes wrist splinting, activity modification, corticosteroid injections, and carpal tunnel release surgery.""",
            metadata={"source": "carpal_tunnel.txt", "topic": "carpal_tunnel", "category": "nerves"}
        ),
        
        # Tendinitis
        Document(
            page_content="""Tendinitis is inflammation of tendons, often caused by repetitive movements or overuse. 
            Common types include tennis elbow, golfer's elbow, and Achilles tendinitis. 
            Treatment includes rest, ice, compression, elevation (RICE), physical therapy, and anti-inflammatory medications. 
            Prevention through proper technique and gradual activity progression is important.""",
            metadata={"source": "tendinitis.txt", "topic": "tendinitis", "category": "soft_tissue"}
        ),
        
        # Bursitis
        Document(
            page_content="""Bursitis is inflammation of bursae, fluid-filled sacs that cushion joints. 
            Common locations include shoulder, elbow, hip, and knee. 
            Causes include repetitive movements, trauma, and inflammatory conditions. 
            Treatment includes rest, ice, anti-inflammatory medications, and corticosteroid injections.""",
            metadata={"source": "bursitis.txt", "topic": "bursitis", "category": "soft_tissue"}
        ),
        
        # Fibromyalgia
        Document(
            page_content="""Fibromyalgia is a chronic condition characterized by widespread musculoskeletal pain, fatigue, and tender points. 
            It's often accompanied by sleep disturbances, mood issues, and cognitive difficulties. 
            Treatment includes medication, exercise, stress management, and cognitive behavioral therapy. 
            A multidisciplinary approach is often most effective.""",
            metadata={"source": "fibromyalgia.txt", "topic": "fibromyalgia", "category": "systemic"}
        ),
        
        # Sports Injuries
        Document(
            page_content="""Sports injuries include sprains, strains, fractures, and overuse injuries. 
            Common injuries include ACL tears, rotator cuff injuries, and stress fractures. 
            Treatment follows the RICE protocol initially, followed by physical therapy and gradual return to activity. 
            Prevention through proper training, equipment, and technique is essential.""",
            metadata={"source": "sports_injuries.txt", "topic": "sports_injuries", "category": "trauma"}
        )
    ]
    
    with st.spinner("Setting up medical knowledge base..."):
        # Split documents into chunks
        text_splitter = RecursiveCharacterTextSplitter(
            chunk_size=1000,
            chunk_overlap=200
        )
        split_docs = text_splitter.split_documents(medical_documents)
        
        # Create vector store
        st.session_state.vectorstore = Chroma.from_documents(
            documents=split_docs,
            embedding=st.session_state.embeddings,
            collection_name="medical_knowledge_base"
        )
    
    return st.session_state.vectorstore

def retrieve_context(query, k=3):
    """Retrieve relevant context from the knowledge base."""
    if st.session_state.vectorstore is None:
        return ""
    
    try:
        docs = st.session_state.vectorstore.similarity_search(query, k=k)
        context = "Relevant information from medical knowledge base:\n\n"
        for i, doc in enumerate(docs, 1):
            context += f"Document {i} (Source: {doc.metadata.get('source', 'Unknown')}):\n"
            context += f"{doc.page_content.strip()}\n\n"
        return context
    except Exception as e:
        st.error(f"Error retrieving context: {e}")
        return ""

def generate_response(question, context, pubmed_articles):
    """Generate a comprehensive medical response."""
    
    # Extract knowledge base info
    kb_info = context if "medical knowledge base" in context else ""
    
    # Create response
    response = f"""## Medical Information for: "{question}"

### 📚 Knowledge Base Information
{kb_info if kb_info else "No specific information found in knowledge base."}

"""
    
    # Add PubMed research if available
    if pubmed_articles:
        response += "### 🔬 Latest Research from PubMed\n"
        for i, article in enumerate(pubmed_articles, 1):
            response += f"""
**Research Article {i}:**
- **Title:** {article['title']}
- **Journal:** {article['journal']}
- **Date:** {article['pub_date']}
- **Authors:** {', '.join(article['authors'][:3])}{' et al.' if len(article['authors']) > 3 else ''}
- **Abstract:** {article['abstract'][:300]}{'...' if len(article['abstract']) > 300 else ''}

"""
    
    # Add general guidance
    response += """
### 💡 General Guidance
This information is for educational purposes only. For specific medical advice, please consult with a qualified healthcare provider.

**When to seek medical attention:**
- Severe or worsening pain
- Pain accompanied by numbness, tingling, or weakness
- Pain that interferes with daily activities
- Pain that doesn't improve with rest and self-care

**Self-care tips:**
- Rest the affected area
- Apply ice for acute injuries, heat for chronic conditions
- Gentle stretching and strengthening exercises
- Maintain good posture and ergonomics
- Consider physical therapy for persistent issues
"""
    
    return response

def main():
    """Main Streamlit application."""
    
    # Header
    st.markdown('<h1 class="main-header">🏥 Medical RAG Assistant</h1>', unsafe_allow_html=True)
    st.markdown('<p style="text-align: center; font-size: 1.2rem; color: #666;">AI-Powered Medical Information with PubMed Research</p>', unsafe_allow_html=True)
    
    # Sidebar
    with st.sidebar:
        st.header("⚙️ Settings")
        
        # PubMed toggle
        include_pubmed = st.checkbox("Include PubMed Research", value=True, help="Search latest medical literature")
        
        # Number of PubMed results
        pubmed_results = st.slider("PubMed Results", 1, 5, 3, help="Number of research articles to retrieve")
        
        # Knowledge base info
        st.header("📚 Knowledge Base")
        st.info("Covers: Back pain, neck pain, arthritis, sciatica, carpal tunnel, tendinitis, bursitis, fibromyalgia, sports injuries, and more.")
        
        # Clear chat button
        if st.button("🗑️ Clear Chat History"):
            st.session_state.chat_history = []
            st.rerun()
    
    # Initialize components
    initialize_embeddings()
    setup_knowledge_base()
    
    # Main chat interface
    st.header("💬 Ask Your Medical Question")
    
    # Chat input
    user_question = st.text_input(
        "What musculoskeletal condition would you like to know about?",
        placeholder="e.g., What causes lower back pain? How to treat sciatica?",
        key="user_input"
    )
    
    # Process question
    if user_question and st.button("🔍 Get Medical Information", type="primary"):
        with st.spinner("Analyzing your question..."):
            
            # Retrieve context
            context = retrieve_context(user_question)
            
            # Get PubMed articles
            pubmed_articles = []
            if include_pubmed:
                pubmed_articles = search_pubmed(user_question, pubmed_results)
            
            # Generate response
            response = generate_response(user_question, context, pubmed_articles)
            
            # Add to chat history
            st.session_state.chat_history.append({
                "question": user_question,
                "response": response,
                "timestamp": time.strftime("%H:%M:%S")
            })
            
            # Clear input
            st.rerun()
    
    # Display chat history
    if st.session_state.chat_history:
        st.header("📋 Chat History")
        
        for i, chat in enumerate(reversed(st.session_state.chat_history)):
            with st.expander(f"Q: {chat['question']} ({chat['timestamp']})", expanded=(i==0)):
                st.markdown(chat['response'])
                
                # Add a delete button for each chat
                col1, col2 = st.columns([1, 10])
                with col1:
                    if st.button("🗑️", key=f"delete_{i}"):
                        st.session_state.chat_history.pop(-(i+1))
                        st.rerun()
    
    # Footer
    st.markdown("---")
    st.markdown("""
    <div style="text-align: center; color: #666; font-size: 0.9rem;">
        <p>⚠️ This tool provides educational information only and should not replace professional medical advice.</p>
        <p>Always consult with a qualified healthcare provider for diagnosis and treatment.</p>
    </div>
    """, unsafe_allow_html=True)

if __name__ == "__main__":
    main() 