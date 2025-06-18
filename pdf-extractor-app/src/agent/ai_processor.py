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
        if not self.api_key:
            raise ValueError("OPENROUTER_API_KEY environment variable is not set")
            
        self.api_url = "https://openrouter.ai/api/v1/chat/completions"
        self.headers = {
            "Authorization": f"Bearer {self.api_key}",
            "HTTP-Referer": "https://github.com/gupadovezi/-pdf-extractor-app-",
            "X-Title": "PDF Extractor App",
            "Content-Type": "application/json"
        }

    def _call_openrouter_api(self, prompt: str) -> str:
        """Call OpenRouter API with the given prompt."""
        try:
            payload = {
                "model": "meta-llama/llama-4-scout:free",
                "messages": [
                    {"role": "system", "content": "You are a helpful AI assistant that extracts and analyzes information from research papers. Always respond with valid JSON."},
                    {"role": "user", "content": prompt}
                ],
                "max_tokens": 4000,
                "temperature": 0.1  # Lower temperature for more consistent JSON output
            }
            
            print(f"Sending request to OpenRouter API...")
            print(f"API URL: {self.api_url}")
            print(f"Headers: {json.dumps({k: v for k, v in self.headers.items() if k != 'Authorization'})}")
            
            response = requests.post(
                self.api_url,
                headers=self.headers,
                json=payload
            )
            
            print(f"Response status code: {response.status_code}")
            
            if response.status_code == 200:
                response_data = response.json()
                if "choices" in response_data and len(response_data["choices"]) > 0:
                    content = response_data["choices"][0]["message"]["content"]
                    print(f"Raw API response: {content[:200]}...")  # Print first 200 chars of response
                    return content
                else:
                    print(f"Unexpected API response format: {json.dumps(response_data)}")
                    raise Exception(f"Unexpected API response format: {response_data}")
            else:
                error_msg = f"Error calling OpenRouter API: {response.status_code}"
                try:
                    error_details = response.json()
                    error_msg += f" - {json.dumps(error_details)}"
                except:
                    error_msg += f" - {response.text}"
                print(f"API error: {error_msg}")
                raise Exception(error_msg)
                
        except Exception as e:
            print(f"API call error: {str(e)}")
            raise Exception(f"Error calling OpenRouter API: {str(e)}")

    def process_text(self, text: str) -> Dict[str, Any]:
        """Process text and extract structured information."""
        try:
            if not text or len(text.strip()) == 0:
                return {
                    "error": "Empty text provided",
                    "raw_text": ""
                }

            print(f"Processing text of length: {len(text)}")
            
            prompt = f"""Extract information from this research paper text and return it as a JSON object with the following structure:
            {{
                "title": "paper title",
                "authors": ["author names"],
                "year": "publication year",
                "journal": "journal name",
                "abstract": "paper abstract",
                "methods": "research methods used",
                "findings": "key findings",
                "conclusions": "main conclusions"
            }}

            Important: Return ONLY the JSON object, no other text or explanation. Do not include markdown formatting or code blocks.

            Paper text:
            {text[:4000]}
            """
            
            response = self._call_openrouter_api(prompt)
            print(f"Received response from API: {response[:200]}...")  # Print first 200 chars of response
            
            # Clean the response to ensure it's valid JSON
            response = response.strip()
            
            # Remove any markdown code block formatting
            if response.startswith("```json"):
                response = response[7:]
            elif response.startswith("```"):
                response = response[3:]
            if response.endswith("```"):
                response = response[:-3]
            response = response.strip()
            
            # Remove any leading/trailing quotes
            response = response.strip('"\'')
            
            try:
                # Try to parse the response as JSON
                data = json.loads(response)
                print(f"Successfully parsed JSON response")
                
                # Validate required fields
                required_fields = ["title", "authors", "year", "journal", "abstract", "methods", "findings", "conclusions"]
                for field in required_fields:
                    if field not in data:
                        print(f"Missing field: {field}")
                        data[field] = ""
                    elif isinstance(data[field], list) and not data[field]:
                        data[field] = []
                    elif isinstance(data[field], str) and not data[field]:
                        data[field] = ""
                
                return data
            except json.JSONDecodeError as e:
                print(f"JSON parsing error: {str(e)}")
                print(f"Raw response: {response}")
                return {
                    "error": f"Failed to parse AI response: {str(e)}",
                    "raw_response": response
                }
                
        except Exception as e:
            print(f"Processing error: {str(e)}")
            return {
                "error": f"Error processing text: {str(e)}",
                "raw_text": text[:200]
            }

    def analyze_findings(self, papers: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze findings across multiple papers."""
        try:
            # Filter out papers with errors
            valid_papers = [p for p in papers if "error" not in p]
            
            if not valid_papers:
                return {
                    "analysis": "No valid papers to analyze",
                    "papers_analyzed": 0
                }
            
            # Create a summary of the papers
            papers_summary = "\n\n".join([
                f"Paper {i+1}:\nTitle: {p.get('title', 'N/A')}\nFindings: {p.get('findings', 'N/A')}"
                for i, p in enumerate(valid_papers)
            ])
            
            prompt = f"""Analyze these research papers and provide insights in JSON format:
            {{
                "common_themes": ["list of common themes"],
                "key_findings": ["list of key findings"],
                "research_gaps": ["list of research gaps"],
                "future_directions": ["list of future research directions"]
            }}

            Papers to analyze:
            {papers_summary}
            
            Important: Return ONLY the JSON object, no other text or explanation. Do not include markdown formatting or code blocks.
            """
            
            analysis = self._call_openrouter_api(prompt)
            
            # Clean and parse the analysis response
            analysis = analysis.strip()
            
            # Remove any markdown code block formatting
            if analysis.startswith("```json"):
                analysis = analysis[7:]
            elif analysis.startswith("```"):
                analysis = analysis[3:]
            if analysis.endswith("```"):
                analysis = analysis[:-3]
            analysis = analysis.strip()
            
            # Remove any leading/trailing quotes
            analysis = analysis.strip('"\'')
            
            try:
                analysis_data = json.loads(analysis)
                return {
                    "analysis": analysis_data,
                    "papers_analyzed": len(valid_papers)
                }
            except json.JSONDecodeError:
                return {
                    "analysis": "Failed to parse analysis response",
                    "raw_analysis": analysis,
                    "papers_analyzed": len(valid_papers)
                }
            
        except Exception as e:
            print(f"Analysis error: {str(e)}")
            return {
                "analysis": f"Error analyzing findings: {str(e)}",
                "papers_analyzed": 0
            } 