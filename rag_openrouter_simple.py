import requests
import json
from langchain_core.documents import Document
from langchain_community.embeddings import HuggingFaceEmbeddings
from langchain_community.vectorstores import Chroma
from langchain.text_splitter import RecursiveCharacterTextSplitter

# Your API configuration
API_KEY = "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"
MODEL_ID = "deepseek/deepseek-v3-base:free"

# Initialize embeddings and vector store
# Using a local embedding model to avoid API key requirements
embeddings = HuggingFaceEmbeddings(model_name="sentence-transformers/all-MiniLM-L6-v2")
vectorstore = None

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

def ask_with_rag(question, temperature=0.7, max_tokens=1024):
    """Ask a question with RAG enhancement."""
    print(f"🤔 Question: {question}")
    
    # Retrieve relevant context
    context = retrieve_context(question)
    
    # Create RAG-enhanced prompt
    if context:
        enhanced_prompt = f"""You are a helpful AI assistant with access to a knowledge base. 
Use the following relevant information to answer the user's question accurately and comprehensively.

{context}

User Question: {question}

Please provide a detailed answer based on the information above. If the information doesn't fully address the question, acknowledge what you know and what might need additional research."""
    else:
        enhanced_prompt = question
    
    # Prepare messages
    messages = [
        {"role": "system", "content": "You are a helpful AI assistant with access to a knowledge base."},
        {"role": "user", "content": enhanced_prompt}
    ]
    
    # Make API call
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
    
    try:
        response = requests.post(
            "https://openrouter.ai/api/v1/chat/completions",
            headers=headers,
            data=json.dumps(payload)
        )
        response.raise_for_status()
        result = response.json()
        
        # Add context information
        if context and "choices" in result:
            result["rag_context"] = context[:500] + "..." if len(context) > 500 else context
        
        return result
    except requests.exceptions.RequestException as e:
        print(f"❌ API call failed: {e}")
        return {"error": str(e)}

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

# Setup knowledge base
setup_knowledge_base(medical_documents, metadata)

# Interactive terminal loop
if __name__ == "__main__":
    print("\n🤖 RAG Medical Assistant Interactive Mode")
    print("Type your question and press Enter. Type 'quit' or 'exit' to stop.\n")
    while True:
        user_input = input("Ask a medical question: ").strip()
        if user_input.lower() in ["quit", "exit"]:
            print("Goodbye!")
            break
        result = ask_with_rag(user_input)
        if "choices" in result:
            print("\nAnswer:", result["choices"][0]["message"]["content"])
            if "rag_context" in result:
                print("\nContext used:", result["rag_context"])
        else:
            print("Sorry, there was an error:", result.get("error")) 