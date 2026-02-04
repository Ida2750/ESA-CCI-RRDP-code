# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 19:21:49 2025

@author: Ida Olsen
"""

import requests
from bs4 import BeautifulSoup
import os

directory = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))) + '/RRDPp/RawData/ASSIST'

print(directory)

os.makedirs(f"{directory}", exist_ok=True)

def get_cruise_csv_links(year):
    base_url = f"https://icewatch.met.no/cruises?year={year}"
    response = requests.get(base_url)
    
    if response.status_code != 200:
        print(f"Failed to fetch data for {year}")
        return []
    
    soup = BeautifulSoup(response.text, 'html.parser')
    cruise_links = [a['href'] for a in soup.select('a') if 'cruise' in a['href']]
    
    csv_links = []
    for cruise in cruise_links:
        cruise_page = requests.get(f"https://icewatch.met.no{cruise}")
        if cruise_page.status_code != 200:
            continue
        
        cruise_soup = BeautifulSoup(cruise_page.text, 'html.parser')
        cruise_name = cruise.split('/')[-1]  # Extract cruise identifier
        for a in cruise_soup.select('a'):
            if 'csv' in a.get('href', '').lower():
                csv_links.append((f"https://icewatch.met.no{a['href']}", cruise_name))
    
    return csv_links

def download_csv(csv_links, year):
    os.makedirs(f"{directory}/{year}", exist_ok=True)
    for link, cruise_name in csv_links:
        filename = f"{link.split('/')[-1][:-4]}-{cruise_name}.csv"  # Ensure unique filename
        response = requests.get(link)
        if response.status_code == 200:
            if 'observations' in filename:
                with open(f"{directory}/{year}/{filename}", 'wb') as file:    
                    file.write(response.content)
                    print(f"Downloaded: {filename}")
        else:
            print(f"Failed to download: {filename}")

# Fetch and download CSVs for 2022 and 2023
years = [2021] #, 2023]
for year in years:
    csv_links = get_cruise_csv_links(year)
    if csv_links:
        download_csv(csv_links, year)
    else:
        print(f"No CSV files found for {year}")