import requests
import json
import os
from typing import List, Dict, Any
from langchain_core.documents import Document
from langchain_openai import OpenAIEmbeddings
from langchain_community.vectorstores import Chroma
from langchain.text_splitter import RecursiveCharacterTextSplitter

class RAGOpenRouter:
    def __init__(self, api_key: str, model_id: str = "deepseek/deepseek-v3-base:free"):
        self.api_key = api_key
        self.model_id = model_id
        self.base_url = "https://openrouter.ai/api/v1/chat/completions"
        self.vectorstore = None
        self.embeddings = OpenAIEmbeddings()
        
    def setup_knowledge_base(self, documents: List[str], metadata: List[Dict] = None):
        """Set up the knowledge base with documents."""
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
        self.vectorstore = Chroma.from_documents(
            documents=split_docs,
            embedding=self.embeddings,
            collection_name="rag_knowledge_base"
        )
        
        print(f"✅ Knowledge base created with {len(split_docs)} document chunks")
        return self.vectorstore
    
    def retrieve_relevant_context(self, query: str, k: int = 3) -> str:
        """Retrieve relevant context from the knowledge base."""
        if not self.vectorstore:
            return ""
        
        try:
            docs = self.vectorstore.similarity_search(query, k=k)
            context = "Relevant information from knowledge base:\n\n"
            for i, doc in enumerate(docs, 1):
                context += f"Document {i} (Source: {doc.metadata.get('source', 'Unknown')}):\n"
                context += f"{doc.page_content.strip()}\n\n"
            return context
        except Exception as e:
            print(f"❌ Error retrieving context: {e}")
            return ""
    
    def create_rag_prompt(self, user_query: str, context: str = "") -> str:
        """Create a RAG-enhanced prompt with context."""
        if context:
            return f"""You are a helpful AI assistant with access to a knowledge base. 
Use the following relevant information to answer the user's question accurately and comprehensively.

{context}

User Question: {user_query}

Please provide a detailed answer based on the information above. If the information doesn't fully address the question, acknowledge what you know and what might need additional research."""
        else:
            return user_query
    
    def call_openrouter_api(self, messages: List[Dict], temperature: float = 0.7, max_tokens: int = 1024) -> Dict:
        """Make API call to OpenRouter."""
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "HTTP-Referer": f"https://openrouter.ai/{self.model_id}",
            "X-Title": "RAG-Enhanced LLM",
            "Content-Type": "application/json"
        }
        
        payload = {
            "model": self.model_id,
            "messages": messages,
            "temperature": temperature,
            "max_tokens": max_tokens,
            "stream": False,
            "top_p": 0.95
        }
        
        try:
            response = requests.post(
                self.base_url,
                headers=headers,
                data=json.dumps(payload)
            )
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"❌ API call failed: {e}")
            return {"error": str(e)}
    
    def ask_with_rag(self, question: str, temperature: float = 0.7, max_tokens: int = 1024) -> Dict:
        """Ask a question with RAG enhancement."""
        print(f"🤔 Question: {question}")
        
        # Retrieve relevant context
        context = self.retrieve_relevant_context(question)
        
        # Create RAG-enhanced prompt
        enhanced_prompt = self.create_rag_prompt(question, context)
        
        # Prepare messages
        messages = [
            {"role": "system", "content": "You are a helpful AI assistant with access to a knowledge base."},
            {"role": "user", "content": enhanced_prompt}
        ]
        
        # Make API call
        result = self.call_openrouter_api(messages, temperature, max_tokens)
        
        # Add context information to result
        if context and "choices" in result:
            result["rag_context"] = context[:500] + "..." if len(context) > 500 else context
        
        return result
    
    def add_documents(self, documents: List[str], metadata: List[Dict] = None):
        """Add new documents to the knowledge base."""
        if not self.vectorstore:
            print("❌ Knowledge base not initialized. Call setup_knowledge_base first.")
            return
        
        # Create Document objects
        docs = []
        for i, content in enumerate(documents):
            meta = metadata[i] if metadata and i < len(metadata) else {"source": f"new_doc_{i}"}
            docs.append(Document(page_content=content, metadata=meta))
        
        # Split documents
        text_splitter = RecursiveCharacterTextSplitter(
            chunk_size=1000,
            chunk_overlap=200
        )
        split_docs = text_splitter.split_documents(docs)
        
        # Add to vector store
        self.vectorstore.add_documents(split_docs)
        print(f"✅ Added {len(split_docs)} new document chunks to knowledge base")


def main():
    """Main function to demonstrate RAG-enhanced OpenRouter usage."""
    
    # Your API key
    API_KEY = "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"
    
    # Initialize RAG system
    rag_system = RAGOpenRouter(API_KEY)
    
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
    rag_system.setup_knowledge_base(medical_documents, metadata)
    
    # Test questions
    test_questions = [
        "How can I treat scoliosis?",
        "What is bracing and when is it used?",
        "Are there surgical options for scoliosis?",
        "What exercises help with scoliosis?",
        "What are alternative treatments for scoliosis?"
    ]
    
    print("\n🧪 Testing RAG-Enhanced System:")
    print("=" * 50)
    
    for i, question in enumerate(test_questions, 1):
        print(f"\nTest {i}: {question}")
        print("-" * 30)
        
        result = rag_system.ask_with_rag(question)
        
        if "choices" in result and result["choices"]:
            response = result["choices"][0]["message"]["content"]
            print(f"🤖 Answer: {response}")
            
            if "rag_context" in result:
                print(f"📚 Context used: {result['rag_context'][:100]}...")
        else:
            print(f"❌ Error: {result.get('error', 'Unknown error')}")
        
        print()
    
    # Interactive mode
    print("\n🤖 Interactive RAG Chat Mode")
    print("=" * 50)
    print("Type 'quit' to exit, 'add' to add new documents")
    
    while True:
        try:
            user_input = input("\n👤 You: ").strip()
            
            if user_input.lower() in ['quit', 'exit', 'bye']:
                print("🤖 Goodbye!")
                break
            elif user_input.lower() == 'add':
                print("📝 Enter new document content (type 'done' when finished):")
                new_docs = []
                while True:
                    doc = input("Document: ").strip()
                    if doc.lower() == 'done':
                        break
                    new_docs.append(doc)
                
                if new_docs:
                    rag_system.add_documents(new_docs)
                    print("✅ Documents added to knowledge base!")
            elif user_input:
                result = rag_system.ask_with_rag(user_input)
                
                if "choices" in result and result["choices"]:
                    response = result["choices"][0]["message"]["content"]
                    print(f"\n🤖 Assistant: {response}")
                    
                    if "rag_context" in result:
                        print(f"\n📚 Context used: {result['rag_context'][:200]}...")
                else:
                    print(f"❌ Error: {result.get('error', 'Unknown error')}")
                    
        except KeyboardInterrupt:
            print("\n\n🤖 Goodbye!")
            break
        except Exception as e:
            print(f"❌ Error: {e}")


if __name__ == "__main__":
    main() 