import os
from typing import List, Dict, Any
from dotenv import load_dotenv
import requests
import json

# Load environment variables
load_dotenv()

class ExtractedData:
    """Schema for extracted data from PDFs"""
    def __init__(self, title="", authors=None, year=None, abstract="", key_findings=None, methodology="", keywords=None):
        self.title = title
        self.authors = authors or []
        self.year = year
        self.abstract = abstract
        self.key_findings = key_findings or []
        self.methodology = methodology
        self.keywords = keywords or []

    def dict(self):
        return {
            "title": self.title,
            "authors": self.authors,
            "year": self.year,
            "abstract": self.abstract,
            "key_findings": self.key_findings,
            "methodology": self.methodology,
            "keywords": self.keywords
        }

class AIProcessor:
    def __init__(self):
        """Initialize the AI processor with OpenRouter API"""
        self.api_key = os.getenv("OPENROUTER_API_KEY")
        self.api_url = "https://openrouter.ai/api/v1/chat/completions"
        self.headers = {
            "Authorization": f"Bearer {self.api_key}",
            "HTTP-Referer": "https://github.com/gupadovezi/-pdf-extractor-app-",  # Required by OpenRouter
            "X-Title": "PDF Extractor App",  # Optional, but good practice
            "Content-Type": "application/json"
        }

    def _call_openrouter_api(self, prompt: str) -> str:
        """Call the OpenRouter API with a prompt"""
        try:
            payload = {
                "model": "meta-llama/llama-4-scout:free",  # Using Llama 4 Scout
                "messages": [
                    {"role": "system", "content": "You are an expert research paper analyzer."},
                    {"role": "user", "content": prompt}
                ],
                "temperature": 0.1,
                "max_tokens": 4000  # Llama 4 Scout has a large context window
            }
            
            response = requests.post(self.api_url, headers=self.headers, json=payload)
            response.raise_for_status()
            
            return response.json()["choices"][0]["message"]["content"]
        except Exception as e:
            print(f"Error calling OpenRouter API: {str(e)}")
            return ""

    def process_text(self, text: str) -> Dict[str, Any]:
        """Process text using AI to extract structured information"""
        try:
            # Create the prompt
            prompt = f"""Extract key information from this research paper text and format it as JSON:
            {text}
            
            Extract the following information:
            - Title
            - Authors (as a list)
            - Year (as a number)
            - Abstract
            - Key findings (as a list)
            - Methodology
            - Keywords (as a list)
            
            Format the response as a valid JSON object with these exact keys."""
            
            # Get response from the model
            response = self._call_openrouter_api(prompt)
            
            # Parse the response
            try:
                data = json.loads(response)
                return ExtractedData(
                    title=data.get("title", ""),
                    authors=data.get("authors", []),
                    year=data.get("year"),
                    abstract=data.get("abstract", ""),
                    key_findings=data.get("key_findings", []),
                    methodology=data.get("methodology", ""),
                    keywords=data.get("keywords", [])
                ).dict()
            except json.JSONDecodeError:
                print("Error parsing JSON response")
                return ExtractedData().dict()
            
        except Exception as e:
            print(f"Error processing text with AI: {str(e)}")
            return ExtractedData().dict()

    def analyze_findings(self, extracted_data: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze multiple papers to find patterns and insights"""
        try:
            # Create a summary prompt
            papers_text = json.dumps(extracted_data, indent=2)
            prompt = f"""Analyze these research papers and provide:
            1. Common themes
            2. Conflicting findings
            3. Research gaps
            4. Future research directions
            
            Papers data:
            {papers_text}
            
            Be objective and data-driven in your analysis."""
            
            # Get analysis from the model
            analysis = self._call_openrouter_api(prompt)
            
            return {
                "analysis": analysis,
                "papers_analyzed": len(extracted_data)
            }
            
        except Exception as e:
            print(f"Error analyzing findings: {str(e)}")
            return {
                "error": str(e),
                "analysis": "",
                "papers_analyzed": 0
            } 