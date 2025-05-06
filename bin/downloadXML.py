#!/usr/bin/env python3

# Import modules
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from webdriver_manager.chrome import ChromeDriverManager
import requests
import re
import os
import time

# Setup Chrome driver
options = webdriver.ChromeOptions()
options.add_argument('--headless')  # Run without opening a browser window
options.add_argument('--no-sandbox')
options.add_argument('--disable-dev-shm-usage')

driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=options)

# Navigate to the HIVDB algorithm updates page
url = "https://hivdb.stanford.edu/page/algorithm-updates/"
driver.get(url)

# Wait a bit for the page to load
time.sleep(3)

# Find all links on the page
links = driver.find_elements(By.TAG_NAME, "a")

# Collect XML links
xml_links = []
for link in links:
    href = link.get_attribute("href")
    if href and href.endswith(".xml"):
        xml_links.append(href)

driver.quit()

if not xml_links:
    print("No XML links found.")
    exit()

# Extract version from filename like HIVdb_2024_03_04.xml
def extract_version(link):
    match = re.search(r'(\d{4}_\d{2}_\d{2})', link)
    return match.group(1) if match else "0000_00_00"

# Sort and get the latest XML link
latest_link = sorted(xml_links, key=extract_version, reverse=True)[0]
xml_filename = os.path.basename(latest_link)

# Download the XML file
response = requests.get(latest_link)
with open(xml_filename, "wb") as f:
    f.write(response.content)