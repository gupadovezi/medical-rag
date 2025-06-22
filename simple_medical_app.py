import streamlit as st
import requests
from datetime import datetime
import json

# Page configuration
st.set_page_config(
    page_title="Medical AI Assistant",
    page_icon="🏥",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #FF6B6B;
        text-align: center;
        margin-bottom: 2rem;
    }
    .chat-message {
        padding: 1rem;
        border-radius: 0.5rem;
        margin-bottom: 1rem;
        border-left: 4px solid #FF6B6B;
    }
    .user-message {
        background-color: #E3F2FD;
        border-left-color: #2196F3;
    }
    .assistant-message {
        background-color: #F3E5F5;
        border-left-color: #9C27B0;
    }
    .stButton > button {
        background-color: #FF6B6B;
        color: white;
        border-radius: 0.5rem;
        border: none;
        padding: 0.5rem 1rem;
    }
    .stButton > button:hover {
        background-color: #FF5252;
    }
</style>
""", unsafe_allow_html=True)

# Initialize session state
if 'chat_history' not in st.session_state:
    st.session_state.chat_history = []

# Medical knowledge base
MEDICAL_KNOWLEDGE = {
    "scoliosis": {
        "description": "Scoliosis is a sideways curvature of the spine that occurs most often during the growth spurt just before puberty.",
        "symptoms": ["Uneven shoulders", "One shoulder blade more prominent", "Uneven waist", "One hip higher than the other"],
        "treatment": ["Observation for mild cases", "Bracing for moderate cases", "Surgery for severe cases"],
        "facts": "Scoliosis affects about 2-3% of the population and is more common in girls than boys."
    },
    "back pain": {
        "description": "Back pain is one of the most common medical problems, affecting 8 out of 10 people at some point.",
        "symptoms": ["Dull aching pain", "Sharp or stabbing pain", "Pain that radiates down the leg", "Muscle stiffness"],
        "treatment": ["Rest and activity modification", "Physical therapy", "Pain medication", "Surgery in severe cases"],
        "facts": "Most back pain is mechanical and improves with conservative treatment within a few weeks."
    },
    "neck pain": {
        "description": "Neck pain can be caused by muscle strain, poor posture, arthritis, or injury.",
        "symptoms": ["Stiffness", "Sharp pain", "Difficulty moving the head", "Headaches"],
        "treatment": ["Rest and gentle stretching", "Physical therapy", "Pain medication", "Posture correction"],
        "facts": "Neck pain is often related to poor posture, especially from looking down at phones or computers."
    },
    "shoulder pain": {
        "description": "Shoulder pain can result from injury, overuse, arthritis, or referred pain from other areas.",
        "symptoms": ["Pain with movement", "Weakness", "Stiffness", "Clicking or popping sounds"],
        "treatment": ["Rest and ice", "Physical therapy", "Anti-inflammatory medication", "Surgery if needed"],
        "facts": "The shoulder is the most mobile joint in the body, making it prone to injury and instability."
    },
    "knee pain": {
        "description": "Knee pain is common and can be caused by injury, arthritis, overuse, or age-related wear.",
        "symptoms": ["Swelling", "Stiffness", "Weakness", "Difficulty walking"],
        "treatment": ["RICE protocol (Rest, Ice, Compression, Elevation)", "Physical therapy", "Pain medication", "Surgery if needed"],
        "facts": "The knee is the largest joint in the body and bears most of our body weight."
    },
    "hip pain": {
        "description": "Hip pain can be caused by arthritis, bursitis, muscle strain, or injury.",
        "symptoms": ["Pain in hip joint", "Pain in groin", "Pain in thigh", "Pain in buttocks"],
        "treatment": ["Rest", "Physical therapy", "Medication", "Hip replacement surgery in severe cases"],
        "facts": "Hip pain often worsens with activity and can significantly impact mobility and quality of life."
    },
    "arthritis": {
        "description": "Arthritis is inflammation of one or more joints, causing pain and stiffness.",
        "symptoms": ["Joint pain", "Stiffness", "Swelling", "Reduced range of motion"],
        "treatment": ["Medication", "Physical therapy", "Lifestyle changes", "Surgery in severe cases"],
        "facts": "There are over 100 types of arthritis, with osteoarthritis and rheumatoid arthritis being the most common."
    },
    "fibromyalgia": {
        "description": "Fibromyalgia is a disorder characterized by widespread musculoskeletal pain.",
        "symptoms": ["Widespread pain", "Fatigue", "Sleep problems", "Memory issues", "Mood changes"],
        "treatment": ["Medication", "Exercise", "Stress management", "Cognitive behavioral therapy"],
        "facts": "Fibromyalgia affects about 2-4% of the population and is more common in women than men."
    },
    "tendonitis": {
        "description": "Tendonitis is inflammation of a tendon, usually caused by overuse or injury.",
        "symptoms": ["Pain near joints", "Stiffness", "Swelling", "Weakness"],
        "treatment": ["Rest", "Ice", "Anti-inflammatory medication", "Physical therapy"],
        "facts": "Tendonitis commonly affects shoulders, elbows, wrists, knees, and ankles."
    },
    "bursitis": {
        "description": "Bursitis is inflammation of the fluid-filled sacs that cushion joints.",
        "symptoms": ["Joint pain", "Swelling", "Stiffness", "Warmth around joint"],
        "treatment": ["Rest", "Ice", "Anti-inflammatory medication", "Corticosteroid injections"],
        "facts": "Bursitis commonly affects shoulders, elbows, hips, and knees."
    }
}

def search_condition(query):
    """Search for medical conditions based on user query"""
    query_lower = query.lower()
    results = []
    
    for condition, info in MEDICAL_KNOWLEDGE.items():
        if condition in query_lower or any(symptom.lower() in query_lower for symptom in info["symptoms"]):
            results.append((condition, info))
    
    return results

def format_response(condition, info):
    """Format the response for a medical condition"""
    response = f"## {condition.title()}\n\n"
    response += f"**Description:** {info['description']}\n\n"
    
    response += "**Common Symptoms:**\n"
    for symptom in info['symptoms']:
        response += f"- {symptom}\n"
    response += "\n"
    
    response += "**Treatment Options:**\n"
    for treatment in info['treatment']:
        response += f"- {treatment}\n"
    response += "\n"
    
    response += f"**Key Facts:** {info['facts']}\n\n"
    
    response += "⚠️ **Important:** This information is for educational purposes only. Always consult with a healthcare professional for proper diagnosis and treatment."
    
    return response

def get_simple_response(query):
    """Get a simple response based on the query"""
    # Search for conditions
    results = search_condition(query)
    
    if results:
        response = ""
        for condition, info in results:
            response += format_response(condition, info) + "\n\n"
        return response
    else:
        # General response for musculoskeletal health
        return """
## General Musculoskeletal Health Information

**Common Causes of Pain:**
- Poor posture
- Overuse or repetitive movements
- Injury or trauma
- Age-related wear and tear
- Inflammatory conditions

**General Prevention Tips:**
- Maintain good posture
- Exercise regularly
- Stretch before and after activity
- Use proper ergonomics
- Maintain a healthy weight

**When to See a Doctor:**
- Severe pain that doesn't improve
- Pain with swelling or redness
- Pain that interferes with daily activities
- Pain that lasts more than a few weeks
- Pain with fever or other symptoms

⚠️ **Important:** This information is for educational purposes only. Always consult with a healthcare professional for proper diagnosis and treatment.
        """

def main():
    # Header
    st.markdown('<h1 class="main-header">🏥 Medical AI Assistant</h1>', unsafe_allow_html=True)
    
    # Sidebar
    with st.sidebar:
        st.header("⚙️ Settings")
        
        # Clear chat button
        if st.button("🗑️ Clear Chat History"):
            st.session_state.chat_history = []
            st.rerun()
        
        # App info
        st.markdown("---")
        st.markdown("### 📊 App Information")
        st.markdown(f"**Chat Messages:** {len(st.session_state.chat_history)}")
        st.markdown("**Knowledge Base:** ✅ Loaded")
        
        # Help
        st.markdown("---")
        st.markdown("### 💡 Tips")
        st.markdown("""
        - Ask about musculoskeletal conditions
        - Be specific with your symptoms
        - The app provides educational information
        - Always consult healthcare professionals
        """)
        
        # Available conditions
        st.markdown("---")
        st.markdown("### 🏥 Available Topics")
        for condition in MEDICAL_KNOWLEDGE.keys():
            st.markdown(f"- {condition.title()}")
    
    # Main chat interface
    col1, col2 = st.columns([3, 1])
    
    with col1:
        # Chat history
        for message in st.session_state.chat_history:
            if message["role"] == "user":
                st.markdown(f"""
                <div class="chat-message user-message">
                    <strong>You:</strong> {message["content"]}
                </div>
                """, unsafe_allow_html=True)
            else:
                st.markdown(f"""
                <div class="chat-message assistant-message">
                    <strong>Assistant:</strong> {message["content"]}
                </div>
                """, unsafe_allow_html=True)
        
        # Input area
        user_input = st.text_area(
            "Ask about musculoskeletal conditions:",
            placeholder="e.g., What are the symptoms of scoliosis?",
            height=100
        )
        
        # Send button
        if st.button("🚀 Send", type="primary"):
            if user_input.strip():
                # Add user message to history
                st.session_state.chat_history.append({
                    "role": "user",
                    "content": user_input,
                    "timestamp": datetime.now().isoformat()
                })
                
                # Get response
                with st.spinner("🤖 Processing your query..."):
                    response = get_simple_response(user_input)
                
                # Add assistant response to history
                st.session_state.chat_history.append({
                    "role": "assistant",
                    "content": response,
                    "timestamp": datetime.now().isoformat()
                })
                
                st.rerun()
    
    with col2:
        st.markdown("### 🏥 Quick Topics")
        quick_topics = [
            "Scoliosis symptoms",
            "Back pain treatment",
            "Neck pain causes",
            "Shoulder pain exercises",
            "Knee pain relief",
            "Arthritis management"
        ]
        
        for topic in quick_topics:
            if st.button(topic, key=f"topic_{topic}"):
                st.session_state.chat_history.append({
                    "role": "user",
                    "content": topic,
                    "timestamp": datetime.now().isoformat()
                })
                
                with st.spinner("🤖 Processing..."):
                    response = get_simple_response(topic)
                
                st.session_state.chat_history.append({
                    "role": "assistant",
                    "content": response,
                    "timestamp": datetime.now().isoformat()
                })
                
                st.rerun()
    
    # Footer
    st.markdown("---")
    st.markdown("""
    <div style='text-align: center; color: #666;'>
        <p>⚠️ <strong>Disclaimer:</strong> This app provides educational information only. 
        Always consult healthcare professionals for medical advice.</p>
        <p>Built with ❤️ using Streamlit and medical knowledge</p>
    </div>
    """, unsafe_allow_html=True)

if __name__ == "__main__":
    main() 