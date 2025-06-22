# 🚀 Medical RAG App Deployment Guide

## Option 1: Streamlit Cloud (Recommended - Free)

### Step 1: Prepare Your Repository
1. **Create a GitHub repository** (if you don't have one)
2. **Push your code** to GitHub:
   ```bash
   git init
   git add .
   git commit -m "Initial commit: Medical RAG App"
   git branch -M main
   git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git
   git push -u origin main
   ```

### Step 2: Deploy on Streamlit Cloud
1. **Go to [share.streamlit.io](https://share.streamlit.io)**
2. **Sign in** with your GitHub account
3. **Click "New app"**
4. **Fill in the details:**
   - **Repository**: `YOUR_USERNAME/YOUR_REPO_NAME`
   - **Branch**: `main`
   - **Main file path**: `rag_web_app.py`
5. **Click "Deploy"**

### Step 3: Your App is Live! 🎉
- **URL**: `https://YOUR_APP_NAME.streamlit.app`
- **Automatic updates** when you push to GitHub
- **Free hosting** with generous limits

---

## Option 2: Railway (Alternative - Free Tier)

### Step 1: Create Railway Account
1. **Go to [railway.app](https://railway.app)**
2. **Sign up** with GitHub
3. **Create new project**

### Step 2: Deploy
1. **Connect your GitHub repository**
2. **Add environment variables** (if needed)
3. **Deploy automatically**

---

## Option 3: Heroku (Paid - More Control)

### Step 1: Create Heroku Account
1. **Go to [heroku.com](https://heroku.com)**
2. **Sign up** and install Heroku CLI

### Step 2: Deploy
```bash
# Install Heroku CLI
brew install heroku/brew/heroku

# Login to Heroku
heroku login

# Create Heroku app
heroku create your-medical-rag-app

# Deploy
git push heroku main

# Open your app
heroku open
```

---

## Option 4: Google Cloud Run (Free Tier)

### Step 1: Setup Google Cloud
1. **Create Google Cloud account**
2. **Enable Cloud Run API**
3. **Install Google Cloud CLI**

### Step 2: Deploy
```bash
# Build and deploy
gcloud run deploy medical-rag-app \
  --source . \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated
```

---

## 🛠️ Files Needed for Deployment

### Required Files:
- ✅ `rag_web_app.py` - Main Streamlit app
- ✅ `requirements.txt` - Dependencies
- ✅ `.streamlit/config.toml` - Streamlit config
- ✅ `README.md` - Project documentation

### Optional Files:
- `Procfile` - For Heroku deployment
- `Dockerfile` - For containerized deployment
- `.dockerignore` - Docker ignore file

---

## 🔧 Environment Variables (if needed)

If you need to set environment variables:

### Streamlit Cloud:
- Go to your app settings
- Add environment variables in the "Secrets" section

### Railway/Heroku:
- Set in the dashboard or via CLI

---

## 📊 Performance Tips

### For Better Performance:
1. **Use caching** for expensive operations
2. **Optimize model loading** (load once, reuse)
3. **Use smaller models** for faster responses
4. **Implement lazy loading** for large datasets

### Monitoring:
- **Streamlit Cloud**: Built-in analytics
- **Railway**: Built-in monitoring
- **Heroku**: Add-on monitoring services

---

## 🚨 Troubleshooting

### Common Issues:

1. **Import Errors**:
   - Check `requirements.txt` has all dependencies
   - Ensure all imports are correct

2. **Memory Issues**:
   - Use smaller models
   - Implement caching
   - Optimize data loading

3. **Timeout Issues**:
   - Reduce model complexity
   - Use async operations
   - Implement progress bars

---

## 🎯 Recommended: Streamlit Cloud

**Why Streamlit Cloud is best:**
- ✅ **Completely free**
- ✅ **Easy deployment** (one-click)
- ✅ **Automatic updates** from GitHub
- ✅ **Built-in analytics**
- ✅ **Custom domains** (optional)
- ✅ **Team collaboration**

**Steps to deploy:**
1. Push code to GitHub
2. Go to share.streamlit.io
3. Connect repository
4. Deploy
5. Share your live URL! 🎉

---

## 📞 Support

If you encounter issues:
- **Streamlit Community**: [discuss.streamlit.io](https://discuss.streamlit.io)
- **GitHub Issues**: Create issue in your repository
- **Documentation**: [docs.streamlit.io](https://docs.streamlit.io) 