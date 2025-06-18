import streamlit as st

st.title("Test App")
st.write("If you can see this, the deployment is working!")

# Test imports
try:
    import pandas as pd
    st.write("✅ Pandas imported successfully")
except ImportError:
    st.write("❌ Pandas import failed")

try:
    import PyPDF2
    st.write("✅ PyPDF2 imported successfully")
except ImportError:
    st.write("❌ PyPDF2 import failed")

try:
    import requests
    st.write("✅ Requests imported successfully")
except ImportError:
    st.write("❌ Requests import failed") 